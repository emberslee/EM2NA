import os
import sys
import math
import time
import numpy as np
np.random.seed(42)
import argparse
import tempfile
from scipy.spatial import cKDTree
from pdbio import (
    chain_names, 
    read_pdb_simple,
    read_atom23_na,
    chain_atom23_to_pdb,
    chains_atom23_to_pdb,
    nuc_type_index_to_1,
)
from sample import (
    get_density,
    filter_with_distance,
    safe_norm,
)
from seqio import (
    read_fasta,
    read_secstr,
    is_dna,
    is_rna,
    align,
    remove_bar,
    remove_bar_simple,
    convert_to_rna_seq,
    get_seq_align_among_candidates,
    format_seqs,
)
from grid import (
    get_aa_type,
)
from tsptrace import trace
from utils import parse_map
from cg2aa import full_atom
from opt_seq_assign_helix import get_ss_by_cssr

def distance(x, y):
    return np.sqrt(np.inner(x-y, x-y))

def split_by_distance(coords, d0):
    segs = []
    seg = []
    for i in range(len(coords)):
        seg.append(coords[i])
        if i + 1 < len(coords):
            v = coords[i] - coords[i+1]
            d = np.sqrt(np.inner(v, v))
            if d > d0:
                if len(seg) > 0:
                    segs.append(seg)
                    seg = []
    if len(seg) > 0:
        segs.append(seg)
    return segs

def extrapolate(p0, p1, p2, d=6.):
    v01 = safe_norm(p1 - p0)
    v12 = safe_norm(p2 - p1)
    cos_theta = np.dot(v01, v12) + 0.1*(np.random.rand()-0.5)*2
    if cos_theta > 1.:
        cos_theta = 1.
    if cos_theta < -1.:
        cos_theta = -1.
    theta = np.arccos(cos_theta)
    axis = np.cross(v01, v12)
    axis = safe_norm(axis)
    R = np.array([
        [np.cos(theta) + axis[0]**2 * (1 - np.cos(theta)), axis[0] * axis[1] * (1 - np.cos(theta)) - axis[2] * np.sin(theta), axis[0] * axis[2] * (1 - np.cos(theta)) + axis[1] * np.sin(theta)],
        [axis[1] * axis[0] * (1 - np.cos(theta)) + axis[2] * np.sin(theta), np.cos(theta) + axis[1]**2 * (1 - np.cos(theta)), axis[1] * axis[2] * (1 - np.cos(theta)) - axis[0] * np.sin(theta)],
        [axis[2] * axis[0] * (1 - np.cos(theta)) - axis[1] * np.sin(theta), axis[2] * axis[1] * (1 - np.cos(theta)) + axis[0] * np.sin(theta), np.cos(theta) + axis[2]**2 * (1 - np.cos(theta))]
    ])
    p3 = p2 + (d + (np.random.rand()-0.5)*2 * 0.2)* np.dot(R, v12)
    return p3


def write_points_as_pdb(filename, res_pos, res_types, bfactors=None, ter=False, chain='A', atom=" P  "):
    assert len(res_types) >= len(res_pos)
    f = open(filename, 'w')
    if bfactors is None:
        bfactors = np.ones(len(res_pos))
    for i in range(len(res_pos)):
        n=i+1
        f.write("ATOM  {:5d} {:4s}   {:1s}{:>2s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n".format(n, atom, res_types[i], chain, n, res_pos[i][0], res_pos[i][1], res_pos[i][2], bfactors[i], bfactors[i]))
        if ter:
            f.write("TER\n")

    if not ter:
        f.write("TER\n")


def query_ball_point(tree, coord, r=6., return_sorted=True):
    inds = tree.query_ball_point(coord, r=r, eps=1e-3)
    if isinstance(inds, int):
        inds = np.asarray([inds], dtype=np.int32)
    data = tree.data
    ds = [np.linalg.norm(data[i] - coord) for i in inds]
    if return_sorted:
        # Sort the indices
        sorted_inds = [i for i in range(len(inds))]
        sorted_inds.sort(key=lambda x:ds[x])
        inds = [inds[x] for x in sorted_inds]
        ds = [ds[x] for x in sorted_inds]
    return ds, inds

def check_path_connectivity(paths, dmax=10., dmin=5.):
    # -1 means we should use less paths
    #  1 means we should use more paths

    # 1. check inter-point distance between two consecutive points
    # if distance in chain is > dmax
    # we should use more paths
    for i, path in enumerate(paths):
        for k in range(len(path)-1):
            if distance(path[k], path[k+1]) > dmax:
                return 1

    # 2. check inter-chain distance
    # if distance between chains is shorted than 6.
    # we should use less paths
    for i in range(len(paths)):
        pi = paths[i]
        for k in range(i+1, len(paths)):
            pk = paths[k]
            d0 = distance(pi[0], pk[0])
            d1 = distance(pi[0], pk[-1])
            d2 = distance(pi[-1], pk[0])
            d3 = distance(pi[-1], pk[-1])
            d = min(d0, d1, d2, d3)
            if d < dmin:
                return -1
    # keep unchanged
    return 0


def convert_atom3_to_atom23(atom3_pos):
    L = len(atom3_pos)
    atom23_pos = np.zeros((L, 23, 3), dtype=np.float32)
    atom23_mask = np.zeros((L, 23), dtype=np.int8)

    #print(atom23_pos[:, [0, 5, 12], :].shape)
    #print(atom3_pos.shape)

    atom23_pos[:, [0, 5, 12], :] = atom3_pos
    atom23_mask[:, [0, 5, 12]] = True
    return atom23_pos, atom23_mask

def fix_atom_pos(atom3, sel_idx=0, dbond=3.5):
    assert sel_idx in [0, 2]
    tree = cKDTree(atom3[:, sel_idx])
    fixed_idxs = set()
    new_atom3 = atom3.copy()
    for i in range(len(atom3)):
        if i in fixed_idxs:
            continue
        idxs = tree.query_ball_point(atom3[i, sel_idx], r=0.50, eps=1e-3)
        idxs = [idx for idx in idxs if idx != i]
        if len(idxs) == 0:
            continue
        fixed_idxs.add(i)
        for idx in idxs:
            if idx in fixed_idxs:
                continue
            curr_c4 = atom3[i, 1]
            clashed_c4 = atom3[idx, 1]
            shift_c4 = clashed_c4 - curr_c4
            updated_pos = atom3[i, sel_idx] + shift_c4
            v = safe_norm(updated_pos - clashed_c4)
            fixed_pos = clashed_c4 + v * dbond
            new_atom3[idx, sel_idx] = fixed_pos
            fixed_idxs.add(idx)
    print("# Fixing {} atoms".format(len(fixed_idxs)))
    return new_atom3

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seq", "-seq", type=str, required=True, help="Input sequences")
    parser.add_argument("--sec", "-sec", type=str, default=None, help="Input secondary structures")
    parser.add_argument("--pmap", "-pmap", type=str, required=True)
    parser.add_argument("--cmap", "-cmap", type=str, required=True)
    parser.add_argument("--nmap", "-nmap", type=str, required=True)
    parser.add_argument("--amap", "-amap", type=str, required=True)
    parser.add_argument("--p", "-p", type=str, required=True)
    parser.add_argument("--c", "-c", type=str, required=True)
    parser.add_argument("--n", "-n", type=str, required=True)
    parser.add_argument("--lkh", "-lkh", type=str, required=True)
    parser.add_argument("--lib", "-lib", type=str, required=True)
    parser.add_argument("--out", "-o", type=str, default="./")
    parser.add_argument("--natype", "-natype", type=str, \
        help="Specify the nucleic-acids type, 'rna/RNA' or 'dna/DNA', if not specified, will be automatically detected")
    # Optional interp
    parser.add_argument("--time_limit", "--time", "-t", type=int, default=172800, help="VRP time limit in seconds")
    parser.add_argument("--pyinterp", action='store_true', help="Interpolation using pyinterp3d")
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    tstart = time.time()

    # Get args
    args = get_args()
    lkh_dir = args.lkh
    lib_dir = args.lib
    lib_bin_dir = lib_dir + '/bin/'
    output_dir = args.out

    # Determine na type
    natype = args.natype
    if natype in ['RNA', 'rna']:
        natype = 'RNA'
    elif natype in ['DNA', 'dna']:
        natype = 'DNA'
    else:
        natype = 'AUTO'

    # Read sequences
    fseq = args.seq
    fsec = args.sec

    seqs = []
    if isinstance(fseq, str) and os.path.exists(fseq):
        seqs = read_fasta(fseq)
        # Format input sequence
        seqs = format_seqs(seqs)

    secs = []
    if isinstance(fsec, str) and os.path.exists(fsec):
        secs = read_secstr(fsec)
 
    # Output
    for i, seq in enumerate(seqs):
        print(f"# Input seq {i} {seq}")
    for i, sec in enumerate(secs):
        print(f"# Input sec {i} {sec}")

    # Filter RNA/DNA
    is_dna_seq = []
    seqs0 = []
    for seq in seqs:
        if is_dna(seq):
            is_dna_seq.append(True)
        else:
            is_dna_seq.append(False)
        # No matter what sequence is input
        # Replace 'T' with 'U'
        seqs0.append(convert_to_rna_seq(seq))
    seqs = seqs0

    # See if all sequences are RNAs or DNAs
    if natype == "AUTO":
        print("# No NA type is input do automatic detection via input sequences")

        #if len(seqs) == 0:
            #natype ="RNA"
            #print("# Detect all sequences to be 'RNA'")
            #pass
        if len(is_dna_seq) > 0:
            if np.all(is_dna_seq):
                natype = "DNA"
                print("# Detect all sequences to be 'DNA'")
            elif np.all( np.logical_not(is_dna_seq) ):
                natype = "RNA"
                print("# Detect all sequences to be 'RNA'")

        if natype == "AUTO":
            # If still AUTO
            print("# Detect all sequences to be both 'DNA' and 'RNA'")

    else:
        print("# Input NA type is '{}' enforce modeling of '{}'".format(natype, natype))
        


    print("# Replace all 'T' with 'U'")
    for i, seq in enumerate(seqs):
        print(f"# Input seq {i} {seq}")
        print(f"# Input seq {i} is dna {is_dna_seq[i]}")

    # Read coords
    fp = args.p
    fc = args.c
    fn = args.n
    pcoords = read_pdb_simple(fp)
    ccoords = read_pdb_simple(fc)
    ncoords = read_pdb_simple(fn)

    # Read map
    fpmap = args.pmap
    fcmap = args.cmap
    fnmap = args.nmap
    famap = args.amap

    if args.pyinterp:
        interp = "py"
    else:
        interp = "f90"
    pmap, porigin, pnxyz, pvsize = parse_map(fpmap, False, None, interp=interp)
    cmap, corigin, cnxyz, cvsize = parse_map(fcmap, False, None, interp=interp)
    nmap, norigin, nnxyz, nvsize = parse_map(fnmap, False, None, interp=interp)
    amap, aorigin, anxyz, avsize = parse_map(famap, False, None, interp=interp)

    pshift = np.copy(porigin)
    cshift = np.copy(corigin)
    nshift = np.copy(norigin)

    # Recluster using radius
    pdens = np.asarray([get_density(pmap, (pcoords[i] - pshift) / pvsize) for i in range(len(pcoords))])
    cdens = np.asarray([get_density(cmap, (ccoords[i] - cshift) / cvsize) for i in range(len(ccoords))])
    ndens = np.asarray([get_density(nmap, (ncoords[i] - nshift) / nvsize) for i in range(len(ncoords))])

    print("# Before filter distance len = {}".format(len(pcoords)))
    print("# Before filter distance len = {}".format(len(ccoords)))
    print("# Before filter distance len = {}".format(len(ncoords)))
    pcoords = filter_with_distance(pcoords, pdens, r0=3.)
    ccoords = filter_with_distance(ccoords, cdens, r0=3.1)
    ncoords = filter_with_distance(ncoords, ndens, r0=3.)
    print("# After  filter distance len = {}".format(len(pcoords)))
    print("# After  filter distance len = {}".format(len(ccoords)))
    print("# After  filter distance len = {}".format(len(ncoords)))
    '''
    for i in range(len(ccoords)):
        for k in range(i+1, len(ccoords)):
            d = distance(ccoords[i], ccoords[k])
            if d < 3.:
                print(i, k, d)
    '''

    # Recalculate
    pinds = np.arange(0, len(pcoords)).tolist()
    cinds = np.arange(0, len(ccoords)).tolist()
    ninds = np.arange(0, len(ncoords)).tolist()

    pdens = np.asarray([get_density(pmap, (pcoords[i] - pshift) / pvsize) for i in range(len(pcoords))])
    cdens = np.asarray([get_density(cmap, (ccoords[i] - cshift) / cvsize) for i in range(len(ccoords))])
    ndens = np.asarray([get_density(nmap, (ncoords[i] - nshift) / nvsize) for i in range(len(ncoords))])

    # Sort
    pinds.sort(key=lambda x : pdens[x], reverse=True)
    cinds.sort(key=lambda x : cdens[x], reverse=True)
    ninds.sort(key=lambda x : ndens[x], reverse=True)

    pcoords = pcoords[pinds]
    ccoords = ccoords[cinds]
    ncoords = ncoords[ninds]

    pdens = pdens[pinds]
    cdens = cdens[cinds]
    ndens = ndens[ninds]

    ptree = cKDTree(pcoords)
    ntree = cKDTree(ncoords)

    # Default
    n_point_per_path = 50 + 40 * math.floor(len(ccoords) / 1000)
    # 1-1000: 50
    # 1000-2000: 90
    # 2000-3000: 130
    

    n_seqs = len(seqs)
    n_search_iter_max = 5
    # For large ribosomes, no need to search a very suitable n_path, so we just set n_search_iter_max to 1
    if 600 <= len(ccoords) < 1000:
        n_search_iter_max = 3
    elif 1000 <= len(ccoords):
        n_search_iter_max = 1

    n_path_min = min(n_seqs, math.floor( len(ccoords) / n_point_per_path ))
    n_path_max = max(n_seqs, math.floor( len(ccoords) / n_point_per_path ))
    n_path_min = max(n_path_min, 1)
    n_path_max = max(n_path_max, 1)
    n_path = (n_path_max + n_path_min) // 2
    print("# Minimum  path is {}".format(n_path_min))
    print("# Maximum  path is {}".format(n_path_max))
    print("# Starting path is {}".format(n_path))

    # Thread using P and C4'
    pccoords = np.concatenate([pcoords, ccoords], axis=0)
    pcdens = np.concatenate([pdens, cdens], axis=0)

    # Start iteratively search
    print("# Start iteratively search", flush=True)

    n_search_iter = 0
    n_seg_min = 3
    n_path_last = None
    csegsx = None
    time_limit = args.time_limit
    print("# VRP Time limit = {:.4f}".format(time_limit))
    while n_search_iter < n_search_iter_max and \
        n_path < n_path_max and \
        n_path > n_path_min and \
        n_path != n_path_last:

        n_path_last = n_path

        print("# Iteration {}/{} thread to {} paths".format(n_search_iter+1, n_search_iter_max, n_path), flush=True)
        n_search_iter += 1
 
        # Trace
        inds0 = trace(pccoords, n_path, lkh_dir, time_limit=time_limit)
        inds0.sort(key=lambda x:len(x), reverse=True)

        segs = []
        psegs = []
        csegs = []

        # After threading to traces (maybe threading by both P and C4')
        # Determine the C4 path by split the trace into P path and C4' path
        for i in range(n_path):
            inds = inds0[i]
            seg = pccoords[inds]
            if len(seg) < 1:
                continue
            segs.append(seg)
            pdens0 = np.asarray([get_density(pmap, (seg[k] - pshift) / pvsize) for k in range(len(seg))])
            cdens0 = np.asarray([get_density(cmap, (seg[k] - cshift) / cvsize) for k in range(len(seg))])
            pseg = []
            cseg = []
            for k in range(len(seg)):
                if pdens0[k] > cdens0[k]:
                    pseg.append(seg[k])
                else:
                    cseg.append(seg[k])
            if len(pseg) > 0:
                psegs.append(pseg)
            if len(cseg) > 0:
                csegs.append(cseg)

        # Check the C4' path
        flag = check_path_connectivity(csegs, dmax=12., dmin=8.)

        # Save the current path
        csegsx = csegs

        # If end here
        if flag == 1:
            n_path = (n_path + n_path_max) // 2
        elif flag == -1:
            n_path = (n_path + n_path_min) // 2
        else:
            print("# Current paths are good, end search")
            break

    # Final check on csegsx
    if csegsx is None:
        inds0 = trace(pccoords, n_path_max, lkh_dir, time_limit=time_limit)
        inds0.sort(key=lambda x:len(x), reverse=True)
        segs = []
        psegs = []
        csegs = []
        # After threading to traces (maybe threading by both P and C4')
        # Determine the C4 path by split the trace into P path and C4' path
        for i in range(n_path_max):
            inds = inds0[i]
            seg = pccoords[inds]
            if len(seg) < 1:
                continue
            segs.append(seg)
            pdens0 = np.asarray([get_density(pmap, (seg[k] - pshift) / pvsize) for k in range(len(seg))])
            cdens0 = np.asarray([get_density(cmap, (seg[k] - cshift) / cvsize) for k in range(len(seg))])
            pseg = []
            cseg = []
            for k in range(len(seg)):
                if pdens0[k] > cdens0[k]:
                    pseg.append(seg[k])
                else:
                    cseg.append(seg[k])
            if len(pseg) > 0:
                psegs.append(pseg)
            if len(cseg) > 0:
                csegs.append(cseg)
        csegsx = csegs 

    # Final
    print("# Finally thread to {} paths".format(len(csegsx)))
    csegs = csegsx

    # Break at too large di->di+1 for each path
    d0 = 10.
    csegs0 = []
    for i, cseg in enumerate(csegs):
        for cseg0 in split_by_distance(cseg, d0):
            if len(cseg0) < n_seg_min:
                continue
            csegs0.append(cseg0)
    csegs = csegs0  


    # Build P-C4'-N using threaded paths
    # 1. Build P
    dcp = 3.8
    dcn = 3.4

    ext_psegs = []
    for i, cseg in enumerate(csegs):
        # Query n+1 P points to determine the direction
        # Extend one more point at both terminus
        p0, p1, p2 = cseg[2], cseg[1], cseg[0]
        v0 = extrapolate(p0, p1, p2)

        p0, p1, p2 = cseg[-3], cseg[-2], cseg[-1]
        v1 = extrapolate(p0, p1, p2)
        ext_cseg = [v0] + cseg + [v1]
        ext_pseg = []

        # Query nearest p position
        for k in range(len(ext_cseg)-1):
            mid = (ext_cseg[k] + ext_cseg[k+1]) / 2.
            r0 = distance(ext_cseg[k], ext_cseg[k+1]) / 2.
            inds = ptree.query_ball_point(mid, r=5, eps=1e-3)
            # Sort according to distance to r0
            ds = [np.abs(r0 - np.linalg.norm(mid - pcoords[ind])) for ind in inds]
            sorted_inds = [x for _, x in sorted(zip(ds, inds))]
            inds = sorted_inds

            # Get nearest p
            if len(inds) >= 1:
                if 3.5 < np.linalg.norm(ext_cseg[k] - pcoords[inds[0]]) < 4.5:
                    p = pcoords[inds[0]]
                else:
                    v = safe_norm(pcoords[inds[0]] - ext_cseg[k])
                    p = ext_cseg[k] + v * (dcp + 0.2*((np.random.rand()-0.5)*2))
            else:
                v = safe_norm(np.random.rand(3))
                p = ext_cseg[k] + v * (dcp + 0.2*((np.random.rand()-0.5)*2))
            ext_pseg.append(p) 
        ext_pseg = np.asarray(ext_pseg, dtype=np.float32)
        ext_psegs.append(ext_pseg)
        #print(len(ext_pseg))
        #print(len(ext_cseg))


    # 2. Build N with 1-NN query
    ext_nsegs = []
    for i, cseg in enumerate(csegs):
        # Query nearest p position
        ext_nseg = []
        for k in range(len(cseg)):
            d, ind = ntree.query(cseg[k], k=1)
            if 2.5 < d < 4.5:
                n = ncoords[ind]
            else:
                v = safe_norm(ncoords[ind] - cseg[k])
                n = cseg[k] + v * (dcn + 0.2*((np.random.rand())-0.5)*2)
            ext_nseg.append(n)
        ext_nseg = np.asarray(ext_nseg, dtype=np.float32)
        ext_nsegs.append(ext_nseg)
        #print(len(ext_nseg))

    # shapes
    # cseg [N, 3]
    # ext_cseg [N+2, 3]
    # ext_pseg [N+1, 3]
    # ext_nseg [N, 3]


    # Determine chain direction using SVD-alignment (optinal)
    # Determine chain direction using Full-atom and C3'-P bonds (select)
    do3p = 1.6
    atom23_chains = []
    # 2024-0128 add atom3_chains
    atom3_chains = []

    for i, cseg in enumerate(csegs):
        ext_pseg = ext_psegs[i]
        ext_nseg = ext_nsegs[i]
        assert len(ext_pseg) == len(cseg) + 1
        assert len(ext_nseg) == len(cseg)
 
        # Forward  uses ext_pseg[: -1]
        atom3 = np.concatenate(
            [
                np.asarray(ext_pseg[:-1], dtype=np.float32)[:, None, :],
                np.asarray(cseg, dtype=np.float32)[:, None, :],
                np.asarray(ext_nseg, dtype=np.float32)[:, None, :]
            ],
            axis=1
        )
        fatom3 = atom3.copy()

        fatom23, fatom23_mask = full_atom(atom3, lib_dir=lib_dir, flag='f')
        o3_pos = fatom23[:, 8, :]
        p_pos = fatom23[:, 0, :]
        fds = np.sum(np.abs(np.sqrt(np.sum((o3_pos[:-1] - p_pos[1:]) ** 2, axis=1)) - do3p))
        print("# Forward  bond-length is {}".format(fds))

        # Backward uses ext_pseg[1 : ][::-1]
        atom3 = np.concatenate(
            [
                np.asarray(ext_pseg[1: ][::-1], dtype=np.float32)[:, None, :],
                np.asarray(cseg[::-1], dtype=np.float32)[:, None, :],
                np.asarray(ext_nseg[::-1], dtype=np.float32)[:, None, :]
            ],
            axis=1
        )
        batom3 = atom3.copy()

        batom23, batom23_mask = full_atom(atom3, lib_dir=lib_dir, flag='b')
        o3_pos = batom23[:, 8, :]
        p_pos = batom23[:, 0, :]
        bds = np.sum(np.abs(np.sqrt(np.sum((o3_pos[:-1] - p_pos[1:]) ** 2, axis=1)) - do3p))
        print("# Backward bond-length is {}".format(bds))
        
        if fds < bds:
            atom3 = fatom3
            atom23 = fatom23
            atom23_mask = fatom23_mask
            print("# Detected direction is 'forward'")
        else:
            atom3 = batom3
            atom23 = batom23
            atom23_mask = batom23_mask
            print("# Detected direction is 'backward'")

        # 2024-0806
        # Fix too close atom positions
        atom3 = fix_atom_pos(atom3, sel_idx=0, dbond=3.5)
        atom3 = fix_atom_pos(atom3, sel_idx=2, dbond=3.5)
        atom23, atom23_mask = full_atom(atom3, lib_dir=lib_dir)

        # Append
        atom3_chains.append(atom3)
        atom23_chains.append(atom23)

    # 2024-0128
    atom3_chains.sort(key=lambda x:len(x), reverse=True)
    atom23_chains.sort(key=lambda x:len(x), reverse=True)       


    # 2024-0129
    # If no seq is provided
    # 1. Detect helix conformation type A-form or B-form
    # 2. Determine seq using Greedy
    if not len(seqs) > 0:
        print("# No seq is input, do automatically detection of A-form B-form helix")
        from geo import helix_form

        # Using atom3
        # @2
        chains_seqs = []
        for i, atom3 in enumerate(atom3_chains):
            seq1 = [get_aa_type(amap, atom3[i, 1], origin=aorigin) for i in range(len(atom3))]
            seq1 = [nuc_type_index_to_1[x] for x in seq1]
            seq1 = "".join(seq1)
            chains_seqs.append(seq1)
            #print(len(seq1))

        chain_indices = []
        for i, atom3 in enumerate(atom3_chains):
            chain_indices += [i] * len(atom3)
            #print(len([i]*len(atom3)))

        chain_indices = np.asarray(chain_indices, dtype=np.int32)

        atom3_pos = np.concatenate(atom3_chains, axis=0)
        atom23_pos, atom23_mask = convert_atom3_to_atom23(atom3_pos)

        # @1
        pairs = []
        with tempfile.TemporaryDirectory() as temp_dir:
            # Ignore too large structure because it will be too slow
            if len(atom23_pos) < 1000:
                pairs = get_ss_by_cssr(atom23_pos, args.lib, temp_dir)

        chains_is_dna = [False] * len(atom3_chains)

        # Final reconstruct
        chains_atom23_pos = []
        chains_atom23_mask = []
        chains_res_types = []

        # Convert forms to chain flags
        dna_flag = np.zeros(len(atom23_pos), dtype=np.int8)
        n_chains = np.max(chain_indices) + 1

        if natype in ['DNA', 'dna']:
            chains_is_dna = [True]*len(atom3_chains)
            dna_flag = np.ones(len(atom23_pos), dtype=np.int8)
            print("# Setting all chains to be DNA")
        elif natype in ['RNA', 'rna']:
            chains_is_dna = [False]*len(atom3_chains)
            print("# Setting all chains to be RNA")
            dna_flag = np.zeros(len(atom23_pos), dtype=np.int8)
        else:
            # Detection
            forms = helix_form(atom23_pos[:, [0, 5, 12], :], pairs, lib_dir)
            forms = np.asarray(forms, dtype=np.int32)
            #print(forms)

            for i, p in enumerate(pairs):
                if forms[i] == 1:
                    # is DNA
                    dna_flag[p[0]] = 1
                    dna_flag[p[1]] = 1
            #print(dna_flag)



        for i in range(n_chains):
            flag = dna_flag[chain_indices == i]
            #print(flag)

            dna_ratio = np.sum(flag) / (len(flag) + 1e-3)

            is_dna = dna_ratio >= 0.50

            print("# Current chain {} DNA ratio = {:.2f} is DNA = {}".format(i, dna_ratio, is_dna))
            # DNA
            
            chains_is_dna[i] = is_dna
            seq = chains_seqs[i]

            #print(len(atom23_pos[chain_indices == i][:, [0, 5, 12], :]))

            #print(len(seq))

            print(seq)

            ratom23_pos, ratom23_mask = full_atom(
                atom23_pos[chain_indices == i][:, [0, 5, 12], :],
                res_types=[x for x in seq], 
                lib_dir=lib_dir, 
                flag="", 
                natype='DNA' if is_dna else 'RNA',
            )
            #print('dna' if is_dna else 'rna')

            chains_atom23_pos.append(ratom23_pos)
            chains_atom23_mask.append(ratom23_mask)
            chains_res_types.append([x for x in seq])

        #pdb_output_dir = os.path.join(output_dir, "denovo_without_seq.pdb")
        pdb_output_dir = os.path.join(output_dir, "denovo_without_seq.cif")
        chains_atom23_to_pdb(
            filename=pdb_output_dir,
            chains_atom23_pos=chains_atom23_pos,
            chains_atom23_mask=chains_atom23_mask,
            chains_res_types=chains_res_types,
            chains_occupancy=None,
            chains_bfactors=None,
            chains_is_dna=chains_is_dna,
        )
        print("# Model is saved to {}".format(pdb_output_dir))
        
        exit(0)




    # If seq is provided
    # Sequence assignments using given seqs and secs
    n_extra = 3
    stride = 5
    chains_atom23_pos = []
    chains_atom23_mask = []
    chains_res_types = []
    chains_is_dna = []
    n_iter_align = 1
    print("# Iterative aligment n = {}".format(n_iter_align))
    for i, atom23 in enumerate(atom23_chains):

        # Get sequence from map using C4'
        seq1 = [get_aa_type(amap, atom23[i, 5], origin=aorigin) for i in range(len(atom23))]
        seq1 = [nuc_type_index_to_1[x] for x in seq1]
        seq1 = "".join(seq1)

        if len(seqs) > 0:
            # Get seq align
            for i_iter_align in range(n_iter_align):
                seqA, seqB, score, k = get_seq_align_among_candidates(seq1, seqs)
                seq1 = seqA
                print("# Iter {} detected seq is {}".format(i_iter_align, seqA))
                print("# Iter {} best fit seq is {}".format(i_iter_align, seqB))
    
            # Rebuild chain with input natypes
            if natype in ["DNA", "dna"]:
                r_atom23, r_atom23_mask = full_atom(atom23[:, [0, 5, 12], :], seqB, lib_dir=lib_dir, flag="", natype="DNA")
                chain_is_dna = True
            elif natype in ["RNA", "rna"]:
                r_atom23, r_atom23_mask = full_atom(atom23[:, [0, 5, 12], :], seqB, lib_dir=lib_dir, flag="", natype="RNA")
                chain_is_dna = False
            else:
                # natype == "AUTO"
                # Automated detection according to best matched seqs[k]
                if is_dna_seq[k]:
                    r_atom23, r_atom23_mask = full_atom(atom23[:, [0, 5, 12], :], seqB, lib_dir=lib_dir, flag="", natype="DNA")
                    chain_is_dna = True
                else:
                    r_atom23, r_atom23_mask = full_atom(atom23[:, [0, 5, 12], :], seqB, lib_dir=lib_dir, flag="", natype="RNA")
                    chain_is_dna = False
     
            if chain_is_dna:
                print("# Detect current chain to be 'DNA'", flush=True)
            else:
                print("# Detect current chain to be 'RNA'", flush=True)
        else:
            # Use predicted sequence
            seqB = seq1
            if natype in ["RNA", "rna"]:
                r_atom23, r_atom23_mask = full_atom(atom23[:, [0, 5, 12], :], seqB, lib_dir=lib_dir, flag="", natype="RNA")
                chain_is_dna = False
                print("# Detected seq is {}".format(seq1))
                print("# Best fit seq is {}".format(seq1))
                print("# Detect current chain to be 'RNA'", flush=True)
            else:
                r_atom23, r_atom23_mask = full_atom(atom23[:, [0, 5, 12], :], seqB, lib_dir=lib_dir, flag="", natype="DNA")
                chain_is_dna = True
                print("# Detected seq is {}".format(seq1))
                print("# Best fit seq is {}".format(seq1))
                print("# Detect current chain to be 'DNA'", flush=True)

 
        # Write for debug
        #chain_atom23_to_pdb(os.path.join(output_dir, f"nuc_chain_{i}.pdb"), r_atom23, r_atom23_mask, seqB, chain=chain_names[i], is_dna=chain_is_dna)
        chain_atom23_to_pdb(os.path.join(output_dir, f"nuc_chain_{i}.cif"), r_atom23, r_atom23_mask, seqB, chain=chain_names[i], is_dna=chain_is_dna)

        # Append
        chains_atom23_pos.append(r_atom23)
        chains_atom23_mask.append(r_atom23_mask)
        chains_res_types.append(seqB)
        chains_is_dna.append(chain_is_dna)

    
    # 2024-0129
    #pdb_output_dir = os.path.join(output_dir, "raw.pdb")
    pdb_output_dir = os.path.join(output_dir, "raw.cif")
    # Convert atom3 to atom23

    # 2024-0128 output raw
    raw_chains_atom23_pos = []
    raw_chains_atom23_mask = []
    raw_chains_res_types = []
    raw_chains_is_dna = []
    for atom3_pos in atom3_chains:
        atom23_pos, atom23_mask = convert_atom3_to_atom23(atom3_pos)
        raw_chains_atom23_pos.append(atom23_pos)
        raw_chains_atom23_mask.append(atom23_mask)
        raw_chains_res_types.append(['U'] * len(atom3_pos))
        raw_chains_is_dna.append(False)

    chains_atom23_to_pdb(
        pdb_output_dir,
        raw_chains_atom23_pos,
        raw_chains_atom23_mask,
        raw_chains_res_types,
        chains_is_dna=raw_chains_is_dna,
    )
    print("# Model is saved to {}".format(pdb_output_dir))
    

    # Output final structure
    #pdb_output_dir = os.path.join(output_dir, "denovo.pdb")
    pdb_output_dir = os.path.join(output_dir, "denovo.cif")
    chains_atom23_to_pdb(
        pdb_output_dir,
        chains_atom23_pos,
        chains_atom23_mask,
        chains_res_types,
        chains_is_dna=chains_is_dna,
    )

    tend = time.time()
    print("# Model is saved to {}".format(pdb_output_dir))
    print("# Finished modeling in {:.4f} seconds".format(tend - tstart))





    #
