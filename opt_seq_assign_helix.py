import os
import tempfile
import subprocess
import numpy as np
from scipy.spatial import KDTree
from typing import Dict, List
from grid import (
    get_grid_value,
    get_aa_type,
)
from pdbio import (
    nuc_type_index_to_1,
    nuc_type_1_to_index,
    chains_atom23_to_pdb,
)
from geo import (
    pairwise_distances,
)
from seqio import (
    read_lines,
    is_dna,
    convert_to_rna_seq,
    ss_to_pairs,
)
from cg2aa import (
    full_atom,
)


"""
    Detect possible SS in the structure
    0. Should convert sequence from a chain to all 'A'/'U' and another chain to all 'U'/'A'
"""


def write_atoms_as_pdb(filename, atom_pos, atom_mask, res_types=None, atom_names=None, ters=None):
    assert len(atom_pos) == len(atom_mask), "{} {}".format(len(atom_pos), len(atom_mask))
    if res_types is not None:
        assert len(atom_pos) == len(atom_mask), "{} {}".format(len(atom_pos), len(res_types))
    if atom_names is not None:
        assert len(atom_pos) == len(atom_names), "{} {}".format(len(atom_pos), len(atom_names))
    if ters is not None:
        assert len(atom_pos) == len(ters), "{} {}".format(len(atom_pos), len(ter))

    f = open(filename, 'w')
    n = 0
    for i in range(len(atom_pos)):
        assert len(atom_pos[i]) == len(atom_mask[i])

        if res_types is None:
            res_type = "U"
        else:
            res_type = res_types[i]

        for k in range(len(atom_pos[i])):
            if atom_mask[i][k] == False:
                continue

            if atom_names is None:
                atom_name = " P  "
            else:
                atom_name = atom_names[i][k]

            f.write("ATOM  {:5d} {:4s}   {:1s}{:>2s}{:4d}    {:8.3f}{:8.3f}{:8.3f}\n".format(n+1, atom_name, res_type, "A", i+1, atom_pos[i][k][0], atom_pos[i][k][1], atom_pos[i][k][2]))
            n += 1

            if ters is not None and ters[i]:
                f.write("TER\n")

    if ters is None:
        f.write("TER\n")

    f.close()
    
def get_ss_by_cssr(atom23_pos, lib_dir=None, out_dir=None):
    if lib_dir is None:
        raise "# Error Please provide lib_dir"

    prefix = os.path.join(lib_dir, "bin")
    program = os.path.join(prefix, "CSSR")

    ppos = atom23_pos[:, 0,  :]
    cpos = atom23_pos[:, 5,  :]
    npos = atom23_pos[:, 12, :]

    pmask = np.ones(len(ppos), dtype=np.int8)
    cmask = np.ones(len(cpos), dtype=np.int8)
    nmask = np.ones(len(npos), dtype=np.int8)

    # Check if two atoms are too close
    # If two atoms are at same position, cssr will be run failed
    r0 = 0.1
    def exclude_atoms(pos, r0=r0):
        mask = np.ones(len(pos), dtype=np.int8)
        tree = KDTree(pos)
        inds = tree.query_ball_point(pos, r=r0, eps=1e-3)
        for i in range(len(inds)):
            if not mask[i]:
                continue

            ind = np.asarray(inds[i], dtype=np.int32)
            ind = ind[ind != i]
            if len(ind) >= 1:
                src = i
                dsts = ind
                print("# {} is close to {} within radius = {}".format(src, dsts, r0))
                for k in ind:
                    mask[k] = False
        return mask

    pmask = exclude_atoms(ppos)
    cmask = exclude_atoms(cpos)
    nmask = exclude_atoms(npos)

    L = len(atom23_pos)
    paired = np.zeros((L, L), dtype=np.int32)

    #####################################################################################
    # Enumerate frag ii and frag kk and send to cssr and check if they can be paired ####
    #####################################################################################

    atom3_pos = np.stack([ppos, cpos, npos], axis=1)

    nmask = np.zeros(len(atom3_pos), dtype=np.int8)
    atom3_mask = np.stack([pmask, cmask, nmask], axis=1)

    length = 12
    stride = 6
    min_length = 10
    n_print = 20
    for i in range(0, len(atom23_pos), stride):
        ibegin = i
        iend = i + length

        iatom3_pos = atom3_pos[ibegin : iend]
        # ignore too short
        if len(iatom3_pos) < min_length:
            continue

        iatom3_mask = atom3_mask[ibegin : iend]
        ires_types = ["A"] * len(iatom3_pos)
        iatom_names = [[" P  ", " C4'", " N9 "] for i in range(len(iatom3_pos))]
        iters = [False] * len(iatom3_pos)
        iters[-1] = True

        iindex_frag_to_all = [-1] * len(iatom3_pos)
        for ii in range(len(iatom3_pos)):
            iindex_frag_to_all[ii] = ibegin + ii

        if i % n_print == 0:
            print("# {} / {}".format(i, len(atom23_pos)))

        for k in range(i + len(iatom3_pos), len(atom23_pos), stride):
            kbegin = k
            kend = k + length

            katom3_pos = atom3_pos[kbegin : kend]
            # ignore too short
            if len(katom3_pos) < min_length:
                continue

            katom3_mask = atom3_mask[kbegin : kend]
            kres_types = ["U"] * len(katom3_pos)
            katom_names = [[" P  ", " C4'", " N1 "] for i in range(len(katom3_pos))]
            kters = [False] * len(katom3_pos)
            kters[-1] = True

            kindex_frag_to_all = [-1] * len(katom3_pos)
            for kk in range(len(katom3_pos)):
                kindex_frag_to_all[kk] = kbegin + kk

            # Write atom3 to PDB
            fin = os.path.join(out_dir, "3p_i_{:0>6d}_{:0>6d}_k_{:0>6d}_{:0>6d}.pdb".format(ibegin, ibegin + len(iatom3_pos), kbegin, kbegin + len(katom3_pos)))
            write_atoms_as_pdb(
                fin,
                atom_pos=np.concatenate([iatom3_pos, katom3_pos], axis=0),
                atom_mask=np.concatenate([iatom3_mask, katom3_mask], axis=0),
                res_types=ires_types + kres_types,
                atom_names=iatom_names + katom_names,
                ters=iters + kters,
            )

            # Parsing using cssr
            fout = os.path.join(out_dir, "3p_i_{:0>6d}_{:0>6d}_k_{:0>6d}_{:0>6d}.cssr".format(ibegin, ibegin + len(iatom3_pos), kbegin, kbegin + len(katom3_pos)))
            cmd = f"{program} {fin} {fout} -o 1 "
            #print("# Running command {}".format(cmd))
            result = subprocess.run(cmd, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if result.returncode != 0:
                print("# Cannot parse ss by cssr")
                print("# STDOUT")
                print(result.stdout)
                print("# STDERR")
                print(result.stderr)

            index_frag_to_all = iindex_frag_to_all + kindex_frag_to_all

            # Read output
            try:
                ss = read_lines(fout)[0].strip()
                pairs = ss_to_pairs(ss)
                for p in pairs:
                    ii = index_frag_to_all[p[0]]
                    kk = index_frag_to_all[p[1]]

                    paired[ii, kk] += 1
                    paired[kk, ii] += 1

            except Exception as e:
                print("# WARN cannot get secondary structure for fragment start at {} {}".format(i, k))
                print(e)
                continue
            
            #end for eack k fragment
        #end for each i fragment
        #break

    # 0131 Update
    # Ignore the bp check, keep all possible bps

    # Get final pairing using paired dmatrix
    # while maximum(paired) > 0 and left_inds is not empty:
    #   find (i, j)* = argmax(paired)
    #   set paired(i, j) = 0
    #   left_inds.remove(i, j)
    # done
    left_inds = set()
    for i in range(L):
        left_inds.add(i)

    final_pairs = []
    matrix = paired.copy()
    n_iter = 0
    n_iter_max = 5000
    while np.max(matrix) > 0 and left_inds and n_iter < n_iter_max:
        i, k = np.unravel_index(np.argmax(matrix), matrix.shape)

        final_pairs.append((i, k))
        # set
        for l in range(L):
            matrix[i, l] = 0
            matrix[l, i] = 0

        for l in range(L):
            matrix[k, l] = 0
            matrix[l, k] = 0

        left_inds.remove(i)
        left_inds.remove(k)
        
        # add up
        n_iter += 1
        
        if n_iter % 100 == 0:
            print("# Round {}/{}".format(n_iter, n_iter_max))

    final_pairs.sort(key=lambda x:x[0])

    # return
    return final_pairs


def get_ss_from_native_structure(atom23_pos, nat_atom23_pos, nat_ss):
    from seqio import ss_to_pairs
    assert len(nat_ss) == len(nat_atom23_pos), "{} {}".format(len(nat_ss), len(nat_atom23_pos))

    cpos = atom23_pos[:, 5]
    nat_cpos = nat_atom23_pos[:, 5]

    neighbor = dict()
    for i in range(len(nat_cpos)):
        neighbor[i] = -1

    used = set()

    tree = KDTree(cpos)
    for i in range(len(nat_cpos)):
        nn = tree.query(nat_cpos[i], k=1, eps=1e-3)
        if nn[0] < 3.0:
            if nn[1] not in used:
                neighbor[i] = nn[1]

    # Assign secstr
    nat_pairs = ss_to_pairs(nat_ss)
    

    pairs = []
    for nat_p in nat_pairs:
        p0 = neighbor[nat_p[0]]
        p1 = neighbor[nat_p[1]]

        if (not p0 == -1) and (not p1 == -1):
            pairs.append((p0, p1))
    pairs.sort(key=lambda x:x[0])


    #print(pairs)
    return pairs


def pairs_to_matrix(n, pairs):
    ret = np.zeros((n, n), dtype=np.int8)
    for p in pairs:
        i = p[0]
        j = p[1]
        ret[i, j] = 1
        ret[j, i] = 1
    return ret

def evaluate_exact(pred_a, true_a):
    tp_map = np.sign(pred_a)*true_a
    tp = tp_map.sum()
    pred_p = np.sign(pred_a).sum()
    true_p = true_a.sum()
    fp = pred_p - tp
    fn = true_p - tp
    recall = tp/(tp+fn)
    precision = tp/(tp+fp)
    if np.isnan(precision):
        precision = 0
    f1_score = 2*tp/(2*tp + fp + fn)
    return precision, recall, f1_score

