import os
import time
import ortools
import tempfile
import numpy as np
from ortools.sat.python import cp_model
from ortools.linear_solver import pywraplp
from scipy.spatial import KDTree
from typing import Dict, List
from dataclasses import dataclass
from Bio import pairwise2
from seqio import align, multi_align
from dp import (
    dynamic_assign_multi,
)
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
    is_dna,
    convert_to_rna_seq,
    remove_intervals,
    get_seq_align_among_candidates,
)
from cg2aa import (
    full_atom,
)
from opt_seq_assign_helix import (
    get_ss_by_cssr,
)

def pjoin(*args):
    return os.path.join(*args)

@dataclass
class PathSeqAssignInfo:
    frag_length: int=20
    frag_stride: int=2
    chain_idx: int=0
    chain_start: int=0
    chain_end: int=0
    seq_full: str=""
    seq_aligned: str=""
    seq_idx: int=0
    seq_start: int=0
    seq_end: int=0
    score: int=0

def has_overlap(range0, range1):
    start = max(range0[0], range1[0])
    end   = min(range0[1], range1[1])
    return end >= start

def size_overlap(range0, range1):
    start = max(range0[0], range1[0])
    end   = min(range0[1], range1[1])
    return max(0, end - start)

def seq_identity(seq0, seq1):
    """
        0. Use head-to-head identity
        1. Use NW/SW alignment
    """
    seq0 = "".join([x for x in seq0 if x != '-'])
    seq1 = "".join([x for x in seq1 if x != '-'])

    #1. NW/SW alignment
    alignment = align(seq0, seq1, mode="global")[0]
    #1.1 Coverage
    matched_region_length = sum(1 for symbol1, symbol2 in zip(alignment[0], alignment[1]) if symbol1 == symbol2)
    cov0 = matched_region_length / len(seq0)
    cov1 = matched_region_length / len(seq1)
    #1.2 Seq.ID.
    identical_symbols = sum(1 for symbol1, symbol2 in zip(alignment[0], alignment[1]) if symbol1 == symbol2)
    ID = identical_symbols / len(alignment[0])
    return ID, cov0, cov1

def log_constraint(assign0, assign1, description=""):
    pass


def add_constraint(path_seq_assign_infos : List[PathSeqAssignInfo], frag_length=20, n_print=50):
    """
        Use aligned positions to build constraint
        constraint[i, j] = 1 means the assignment between i and j is conflicted
        and only i or j can be reserved
    """
    constraint_types = [0] * 8

    L = len(path_seq_assign_infos)
    if L >= 10000:
        n_print = 200
    constraint = np.zeros((L, L)).astype(np.int32)
    for i in range(L):
        i_assign = path_seq_assign_infos[i]
        i_frag_length = i_assign.frag_length

        for k in range(i+1, L):
            k_assign = path_seq_assign_infos[k]
            k_frag_length = k_assign.frag_length

            # 0 If two assign is from the same chain-seq
            if  i_assign.chain_idx == k_assign.chain_idx and \
                i_assign.seq_idx   == k_assign.seq_idx:
                constraint[i, k] = 1
                constraint[k, i] = 1
                constraint_types[5] += 1
                log_constraint(i_assign, k_assign, "type5")
                continue

            # 1. If two fragment in same path
            in_same_path = False
            in_same_path = (i_assign.chain_idx == k_assign.chain_idx)
            has_frag_overlap = False

            if in_same_path:
                s0, e0 = i_assign.chain_start, i_assign.chain_end
                s1, e1 = k_assign.chain_start, k_assign.chain_end
                # [s0, ...., e0] | [s1, ..., e1]
                union_s = min(s0, s1)
                union_e = max(e0, e1)
                if union_s - union_e < (i_frag_length + k_frag_length):
                    has_frag_overlap = True
                # If have overlap, we should add extra constarint
                if has_frag_overlap:
                    # [s0, ...., e0] & [s1, ..., e1]
                    inter_s = max(s0, s1)
                    inter_e = min(e0, e1)
                    # Determine which is first
                    if union_s == s0:
                        assign0, assign1 = i_assign, k_assign
                    else:
                        assign0, assign1 = k_assign, i_assign
                    # 1.1 Should from same sequence
                    if assign0.seq_idx != assign1.seq_idx:
                        constraint[i, k] = 1
                        constraint[k, i] = 1
                        constraint_types[0] += 1
                        log_constraint(assign0, assign1, "type0")
                        continue
                    # 1.2 Should have overlapped seq-assign intervals
                    if not has_overlap(
                            [assign0.seq_start, assign0.seq_end], 
                            [assign1.seq_start, assign1.seq_end],
                        ):
                        constraint[i, k] = 1
                        constraint[k, i] = 1
                        constraint_types[1] += 1
                        log_constraint(assign0, assign1, "type1")
                        continue
                    # Get the aligned seq (without gap) for both fragment
                    align0_begin = inter_s - union_s
                    align0_useful_count = len([c for c in assign0.seq_aligned[:align0_begin] if c != "-"])
                    align0_interval = [assign0.seq_start + align0_useful_count, assign0.seq_end]
                    align0_seq = assign0.seq_aligned[align0_begin:]

                    align1_end = inter_e - inter_s
                    align1_useful_count = len([c for c in assign1.seq_aligned[:align1_end] if c != "-"])
                    align1_interval = [assign1.seq_start, assign1.seq_start + align1_useful_count]
                    align1_seq = assign1.seq_aligned[:align1_end]

                    # 1.3 Should have overlap in continues seq
                    # Do not consider the constraint

                    # 1.4 Should have > 80% identity
                    seqid, cov0, cov1 = seq_identity(
                        align0_seq,
                        align1_seq,
                    )
                    if not (seqid >= 0.80 and cov0 >= 0.80 and cov1 >= 0.80):
                        constraint[i, k] = 1
                        constraint[k, i] = 1
                        constraint_types[3] += 1
                        log_constraint(assign0, assign1, "type3 seq.id={:.2f} seq0={} seq1={}".format(seqid, align0_seq, align1_seq))
                        continue

            # 2. If two fragment NOT in same path
            if not in_same_path:
                assign0, assign1 = i_assign, k_assign

                if assign0.seq_idx != assign1.seq_idx:
                    # 2.0 Different sequence, acceptable
                    pass
                else:
                    # 2.1 Same sequence, should have non-overlapped sequence interval
                    min_frag_length = min(i_frag_length, k_frag_length)
                    max_frag_length = max(i_frag_length, k_frag_length)
                    l = size_overlap(
                            [assign0.chain_start, assign0.chain_end],
                            [assign1.chain_start, assign1.chain_end],
                        )
                    if l >= 0.2 * max_frag_length:
                        constraint[i, k] = 1
                        constraint[k, i] = 1
                        constraint_types[4] += 1
                        continue

        if i % n_print == 0:
            print("# Report adding constraint {} / {}".format(i, L), flush=True)

    #for i in range(6):
    #    print("# Constraint type {} count {}".format(i, constraint_types[i]))

    return constraint


def score_assign(path_seq_assign_infos : List[PathSeqAssignInfo], rescale=10.):
    """
        Use seqeence alignment score 
    """
    L = len(path_seq_assign_infos)
    S = np.zeros(L).astype(np.int32)
    for k in range(L):
        assign = path_seq_assign_infos[k]
        s = int(assign.score * rescale)
        S[k] = s
    return S



def solve_seq_assign(path_seq_assign_infos, constraint):
    """
        The constraint is a [N, N] matrix, N is the no. of fragments
        It defines whether fragment i, j has clashes
    """
    model = cp_model.CpModel()

    x = []
    for i in range(len(constraint)):
        x.append(model.NewBoolVar(f'x[{i}]'))

    # Node constraint
    print("# Calculating node constraints", flush=True)
    for i in range(len(constraint)):
        model.Add(sum(x[j] for j in range(len(constraint)) if constraint[i][j]==1)==0).OnlyEnforceIf(x[i])

    # Score matrix
    print("# Calculating score matrix", flush=True)
    score_array = score_assign(path_seq_assign_infos)
    score_array = score_array.astype(int)

    # Define optimizing target
    print("# Defining optimizing target", flush=True)
    objective_terms = []
    for i in range(len(constraint)):
        objective_terms.append(score_array[i]* x[i]) #Sum of Raw Scores
    model.Maximize(sum(objective_terms))

    # Solve
    print("# Start Solving...", flush=True)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = 1800
    solver.parameters.num_search_workers = 4
    solver.parameters.log_search_progress = False
    solution_printer = cp_model.ObjectiveSolutionPrinter()

    # Print status
    status = solver.SolveWithSolutionCallback(model, solution_printer)
    results = []
    if status == pywraplp.Solver.OPTIMAL or pywraplp.Solver.FEASIBLE:
        for i in range(len(constraint)):
            if solver.BooleanValue(x[i]):
                results.append(i)
    return results


if __name__ == '__main__':
    tstart = time.time()
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", "-p", required=True, type=str, help="Input fragments")
    #parser.add_argument("--raw", "-r", required=True, type=str, help="Input raw fragments")
    parser.add_argument("--aamap", "-a", required=True, type=str, help="Input aa grid")
    parser.add_argument("--aaprob", required=True, type=str, help="Input 4-class aa probability")
    parser.add_argument("--seq", "-s", required=True, type=str, help="Input sequence")
    parser.add_argument("--lib", "-l", required=True, type=str, help="Input directory of lib")
    parser.add_argument("--output", "-o", type=str, default="./", help="Output diretory")
    parser.add_argument("--no_opt_assign", "-no_opt_assign", action='store_true', help="Specify to skip optimizing seq assign")
    parser.add_argument("--no_check_bp", "-no_check_bp", action='store_true', help="Specify to skip checking bp restraint")
    # Optional interp
    parser.add_argument("--pyinterp", action='store_true', help="Interpolation using pyinterp3d")
    args = parser.parse_args()

    lib_dir = args.lib 
    fout = args.output

    fseq = args.seq
    faamap = args.aamap
    fpdb = args.pdb
    #fraw = args.raw

    # Read seq
    from seqio import read_fasta, format_seqs
    seqs = read_fasta(fseq)
    seqs = format_seqs(seqs)
    seqs_is_dna = [False]*len(seqs)
    print("# Found {} input sequence".format(len(seqs)))
    for i in range(len(seqs)):
        print("# {} {}".format(i, seqs[i]))
        if is_dna(seqs[i]):
            seqs_is_dna[i] = True
        seq = seqs[i]
        seq = convert_to_rna_seq(seq)
        seqs[i] = seq

    # 0129
    # If n * DNA duplex set to be all DNA
    if len(seqs_is_dna) % 2 == 0:
        ratio = np.sum(seqs_is_dna) / len(seqs_is_dna)
        print("# DNA seqs ratio is {:.4f}".format(ratio))
        if ratio >= 0.50:
            seqs_is_dna = [True] * len(seqs_is_dna)
            print("# Set to be all DNA")

    # For single sequence
    # Further consider all sequences
    all_same_res = None
    if len(seqs) == 1:
        seq = seqs[0]
        cnts = [0] * 4
        for k in range(len(seq)):
            cnts[nuc_type_1_to_index[seq[k]]] += 1
        ind = np.argmax(cnts)
        ratio = cnts[ind] / np.sum(cnts)
        if ratio >= 0.80:
            all_same_res = nuc_type_index_to_1[ind]
    print("# All same res : {}".format(all_same_res))

    # Read aa map
    from utils import parse_map
    if args.pyinterp:
        interp = "py"
    else:
        interp = "f90"
    aamap0, aaorigin0, aanxyz0, aavsize0 = parse_map(faamap, False, None, interp=interp)

    # Read pdb
    from pdbio import read_atom23_na
    atom23_pos, atom23_mask, res_types, res_indices, chain_indices = read_atom23_na(fpdb)
    print("# Read {} nucleotides".format(len(atom23_pos)), flush=True)

    # Read raw
    """
    raw_atom23_pos, raw_atom23_mask, raw_res_types, raw_res_indices, raw_chain_indices = read_atom23_na(fraw)
    print("# Read {} raw nucleotides".format(len(raw_atom23_pos)), flush=True)
    assert len(atom23_pos) == len(raw_atom23_pos)
    """

    n_chains = np.max(chain_indices) + 1
    chains = []
    for i in range(0, n_chains):
        try:
            chain_atom23_pos = atom23_pos[chain_indices == i]
            chain_atom23_mask = atom23_mask[chain_indices == i]
            if len(chain_atom23_pos) > 0:
                chains.append(chain_atom23_pos)
        except:
            print(f"# WARN cannot find chain {i}")
            continue
    print("# Found {} chains".format(len(chains)))

    # Read aa probability
    faaprob = args.aaprob
    aaprob = np.load(faaprob)
    aaprob = aaprob['prob']

    # Prepare for alignment
    to_frags = False
    FRAG_LENGTH = int(1e5)
    FRAG_STRIDE = int(1e5)
    if to_frags:
        FRAG_LEGNTH = 40
        FRAG_STRIDE = 5

    frags_aaprobs = []
    frags_dmat = []
    frags_cpos = []
    frags_atom23_pos = []
    frags_cpos = []

    for _, atom23_pos in enumerate(chains):
        # atom23
        frag_atom23_pos = atom23_pos
        frags_atom23_pos.append(frag_atom23_pos)

        # cpos
        cpos = frag_atom23_pos[:, 5]
        npos = frag_atom23_pos[:, 12]

        frags_cpos.append(cpos)
    
        # aaprob
        aaprobs = []
        for i in range(len(cpos)):
            probs = []
            for k in range(4):
                prob = get_grid_value(aaprob[k], (npos[i]-aaorigin0)/aavsize0)
                probs.append(prob)
            aaprobs.append(probs)
        frags_aaprobs.append(aaprobs)
    
        # dmat
        dmat = pairwise_distances(cpos, cpos)
        frags_dmat.append(dmat)

    infos = []
    assign_infos = []
    for i, atom23_pos in enumerate(frags_atom23_pos):
        aaprobs = frags_aaprobs[i]
        dmat = frags_dmat[i]

        for k, seq in enumerate(seqs):
            if len(seq) / len(atom23_pos) < 1.00:
                print("# Ignore, try to assign seq of len {} to frag of len {}".format(len(seq), len(atom23_pos)))
                continue

            seq = [nuc_type_1_to_index[x] for x in seq]
            max_score, match_seq, match_interval = dynamic_assign_multi(aaprobs, dmat, seq, save_path=None, top=10)

            for l in range(len(max_score)):
                infos.append(
                    PathSeqAssignInfo(
                        frag_length=len(atom23_pos),
                        frag_stride=int(1e5),

                        chain_idx=i,
                        chain_start=0,
                        chain_end=len(atom23_pos)-1,

                        seq_full=seqs[k],
                        seq_aligned=match_seq[l],

                        seq_idx=k,
                        seq_start=match_interval[l][0],
                        seq_end=match_interval[l][1],
                        score=int(max_score[l]),
                    )
                )

    print("# In total we have {} path-seq assignments".format(len(infos)))

    # Add constraint
    constraint = add_constraint(infos)
    
    # Check if chain is reassigned, if not output anyway
    r_chains_is_reassigned = []
    
    # Solve assign
    result = solve_seq_assign(infos, constraint)
    print("# Result ", result)

    # Check result
    #print("# Check optimize result")
    r_chains_res_types = []
    r_chains_atom23_pos = []
    r_chains_atom23_mask = []
    r_chains_is_dna = []

    r_chains_atom3_pos = []

    # Greedy assign the left chains
    left_chains_res_types = []
    left_chains_atom23_pos = []
    left_chains_atom23_mask = []
    left_chains_is_dna = []

    left_chains_atom3_pos = []

    if not len(result) > 0:
        print("# Cannot optimize the sequence assignment")

    # Merge result
    realign = True
    assigned_seq_ranges = dict()
    for kk in range(len(seqs)):
        assigned_seq_ranges[kk] = list()

    for idx in result:
        info = infos[idx]

        # To avoid repeating (impossible to happen but who knows)
        if info.chain_idx in r_chains_is_reassigned:
            continue
    
        assigned_seq_ranges[info.seq_idx].append([info.seq_start, info.seq_end])

        r_chains_is_reassigned.append(info.chain_idx)
       
        aaprobs = frags_aaprobs[info.chain_idx]
        dmat = frags_dmat[info.chain_idx]
        seq_aligned = info.seq_aligned
        
        if realign:
            seq = seqs[info.seq_idx]
            start = info.seq_start
            end   = info.seq_end
            seq_str = seq[start : end]

            seq1 = np.argmax(np.asarray(aaprobs), axis=1).tolist()
            seq1 = "".join([nuc_type_index_to_1[x] for x in seq1])
            seqA, seqB, score, kk = get_seq_align_among_candidates(seq1, [seq_str], setting=1)

            # post process
            seqBx = ""
            for kk in range(len(seqA)):
                if seqB[kk] == '-':
                    seqBx += seqA[kk]
                else:
                    base = np.argmax(aaprobs[kk])
                    if base != nuc_type_1_to_index[seqB[kk]]:
                        print("# Mismatched {} {}".format(base, nuc_type_1_to_index[seqB[kk]]))
                    if aaprobs[kk][base] >= 0.80:
                        print("# Use {} instead".format(base))
                        seqBx += nuc_type_index_to_1[base]
                    else:
                        seqBx += seqB[kk]
            seqB = seqBx

            # final
            seq_aligned = seqB

 
        r_chains_is_dna.append(seqs_is_dna[info.seq_idx])
        natype = "dna" if seqs_is_dna[info.seq_idx] else "rna"
    
        r_res_types = [x if x != '-' else 'A' for x in seq_aligned]
        r_chains_res_types.append(r_res_types)
    
        # Use new alignment anyway
        atom23_pos = frags_atom23_pos[info.chain_idx]
        atom3_pos = atom23_pos[:, [0, 5, 12], :]
        r_chains_atom3_pos.append(atom3_pos)

        r_atom23_pos, r_atom23_mask = full_atom(atom3_pos, res_types=r_res_types, lib_dir=lib_dir, flag="", natype=natype)
    
        r_chains_atom23_pos.append(r_atom23_pos)
        r_chains_atom23_mask.append(r_atom23_mask)

    # Exclude assigned target seqs
    exclude_assigned_seqs = False
    left_seqs = []
    left_seqs_is_dna = []

    if not exclude_assigned_seqs:
        left_seqs = seqs
        left_seqs_is_dna = seqs_is_dna
        print("# No exclude assigned seqs")
    else:
        pass
        for k in range(len(seqs)):
            left_intervals = remove_intervals(len(seqs[k]), assigned_seq_ranges[k])
            print(left_intervals)
            for interval in left_intervals:
                s, e = interval
                # Ignore too short seq
                if e - s < -1:
                    continue
                left_seqs.append(seqs[k][s : e])
                left_seqs_is_dna.append(seqs_is_dna[k])

        print("# Exclude assigned seqs {} sub-seqs left".format(len(left_seqs)))

    print("# Greedy assign for left chains")
    for i in range(len(frags_atom23_pos)):
        if i in r_chains_is_reassigned:
            continue
        print("# Greedy assign left chain {} left".format(i))
        aaprobs = frags_aaprobs[i]
        dmat = frags_dmat[i]

        l = len(frags_atom23_pos[i])

        seq1 = [get_aa_type(aamap0, frags_atom23_pos[i][k, 5], origin=aaorigin0) for k in range(l)]
        seq1 = [nuc_type_index_to_1[x] for x in seq1]
        seq1 = "".join(seq1)

        seqA, seqB, score, kk = get_seq_align_among_candidates(seq1, left_seqs, setting=1)
        assert len(seqA) == len(seqB)


        best_match_seq = seqB
        best_max_score = score
        best_match_seq_k = kk
        best_match_seq_start = -1
        best_match_seq_end = -1
        best_match_seq_is_dna = left_seqs_is_dna[kk]


        # If still unassigned
        aatypes = np.argmax(np.asarray(aaprobs, dtype=np.float32), axis=-1)
        aatypes = [nuc_type_index_to_1[x] for x in aatypes]

        if best_match_seq is None:
            print("# Use prediction directly from map")
            # Use prediction from map
            best_match_seq = aatypes
            best_match_seq = "".join(best_match_seq)
        else:
            print("# Use best match seq = {} range = [{}, {}) score = {:.2f}".format(best_match_seq_k, best_match_seq_start, best_match_seq_end, best_max_score))

        # Append
        natype = "dna" if best_match_seq_is_dna else "rna"
        best_match_seq_list = [best_match_seq[i] if best_match_seq[i] != '-' else aatypes[i] for i in range(len(best_match_seq))]

        atom23_pos = frags_atom23_pos[i]
        atom3_pos = atom23_pos[:, [0, 5, 12], :]
        left_chains_atom3_pos.append(atom3_pos)

        atom23_pos, atom23_mask = full_atom(atom3_pos, res_types=best_match_seq_list, lib_dir=lib_dir, flag="", natype=natype)

        left_chains_atom23_pos.append(atom23_pos)
        left_chains_atom23_mask.append(atom23_mask)
        left_chains_is_dna.append(best_match_seq_is_dna)
        left_chains_res_types.append(best_match_seq_list)


    # Concat final chains
    final_chains_res_types   = r_chains_res_types   + left_chains_res_types
    final_chains_atom23_pos  = r_chains_atom23_pos  + left_chains_atom23_pos
    final_chains_atom23_mask = r_chains_atom23_mask + left_chains_atom23_mask
    final_chains_is_dna      = r_chains_is_dna      + left_chains_is_dna

    final_chains_atom3_pos   = r_chains_atom3_pos   + left_chains_atom3_pos

    assert len(final_chains_res_types) == len(final_chains_atom23_pos) == \
           len(final_chains_atom23_mask) == len(final_chains_is_dna) == len(final_chains_atom3_pos)
          
    # Sort the chains
    print("# Sorting the chain length in descending order")
    inds = [i for i in range(len(final_chains_atom23_pos))]
    inds.sort(key=lambda x:len(final_chains_atom23_pos[x]), reverse=True)

    final_chains_res_types   = [final_chains_res_types[x] for x in inds]
    final_chains_atom23_pos  = [final_chains_atom23_pos[x] for x in inds]
    final_chains_atom23_mask = [final_chains_atom23_mask[x] for x in inds]
    final_chains_is_dna      = [final_chains_is_dna[x] for x in inds]
    
    final_chains_atom3_pos  = [final_chains_atom3_pos[x] for x in inds]

    # For single sequence
    if all_same_res is not None:
        final_chains_res_types = [[all_same_res] * len(x) for x in final_chains_atom23_pos]

    # Write chains before base-type checking
    #fpdbout = pjoin(fout, "reassign.pdb")
    fpdbout = pjoin(fout, "reassign.cif")
    print("# Output final chains to {}".format(fpdbout))
    chains_atom23_to_pdb(
        fpdbout, 
        final_chains_atom23_pos, 
        final_chains_atom23_mask, 
        final_chains_res_types, 
        chains_is_dna=final_chains_is_dna
    )


    # Check if skip bp check
    if args.no_check_bp:
        print("# No check base pair restraint")
        #fpdbout = pjoin(fout, "reassign_helix.pdb")
        fpdbout = pjoin(fout, "reassign_helix.cif")
        print("# Output final chains to {}".format(fpdbout))
        chains_atom23_to_pdb(
            fpdbout, 
            final_chains_atom23_pos, 
            final_chains_atom23_mask, 
            final_chains_res_types, 
            chains_is_dna=final_chains_is_dna
        )
        tend = time.time()
        print("# Time consuming {:.4f}".format(tend - tstart))
        exit(0)


    #######################################
    # Start check bp and opt helix assign #
    #######################################
    chain_indices = []
    for i, chain_atom23_pos in enumerate(final_chains_atom23_pos):
        chain_indices += [i] * len(chain_atom23_pos)
    chain_indices = np.asarray(chain_indices, dtype=np.int32)

    res_types     = []
    for _, chain_res_types in enumerate(final_chains_res_types):
        res_types += [nuc_type_1_to_index[x] for x in chain_res_types]
    res_types = np.asarray(res_types, dtype=np.int32)

    atom23_pos    = np.concatenate(final_chains_atom23_pos, axis=0)
    atom23_mask   = np.concatenate(final_chains_atom23_mask, axis=0)

    atom3_pos = np.concatenate(final_chains_atom3_pos, axis=0)

    # Recalculate aaprobs
    aaprobs = []
    npos = atom23_pos[:, 12]
    for i in range(len(npos)):
        probs = np.zeros(4)
        div = np.zeros(4)
        ranges = [0]
        for k in range(4):
            x, y, z = npos[i].tolist()

            for xd in ranges:
                for yd in ranges:
                    for zd in ranges:
                        xpos, ypos, zpos = x + xd, y + yd, z + zd
                        pos = np.asarray([xpos, ypos, zpos])
                        prob = get_grid_value(aaprob[k], (pos-aaorigin0)/aavsize0)
                        probs[k] += prob
                        div[k] += 1
        probs = (probs / div).tolist()
        aaprobs.append(probs)
    aaprobs = np.asarray(aaprobs)

    # Detect base pairs
    pairs = []
    with tempfile.TemporaryDirectory() as temp_dir:
        # Ignore too large structure because it will be too slow
        if len(atom23_pos) < 1000:
            pairs = get_ss_by_cssr(atom23_pos, lib_dir, temp_dir)
            for p in pairs:
                print("# Detected bp with index {} - {}".format(p[0], p[1]))


    # Post process
    watson_crick_pair = {
        'A': 'U',
        'U': 'A',
        'C': 'G',
        'G': 'C',
    }
    res_types0 = res_types.copy()
    for p in pairs:
        base0 = res_types[p[0]]
        base1 = res_types[p[1]]

        # 1. use assigned base type
        prob0 = aaprobs[p[0]][base0]
        prob1 = aaprobs[p[1]][base1]
        probs = [prob0, prob1]

        ind = np.argmax(probs)
        rbase0 = base0
        rbase1 = base1
        if probs[ind] >= 0.50:
            if ind == 0:
                rbase1 = nuc_type_1_to_index[watson_crick_pair[nuc_type_index_to_1[base0]]]
            else:
                rbase0 = nuc_type_1_to_index[watson_crick_pair[nuc_type_index_to_1[base1]]]
            print("# Correct base type {} {} to {} {}".format(base0, base1, rbase0, rbase1))

        res_types[p[0]] = rbase0
        res_types[p[1]] = rbase1

    res_types = [nuc_type_index_to_1[x] for x in res_types]

    # For single sequence
    if all_same_res is not None:
        res_types = [all_same_res] * len(res_types)

    # Reconstruct full-atom
    from build import full_atom
    n_chains = np.max(chain_indices) + 1
    chains_atom23_pos = []
    chains_atom23_mask = []
    chains_res_types = []
    chains_is_dna = []
    print("# Num chains = {}".format(n_chains))

    #print(chain_indices)
    #print(raw_chain_indices)
   
    for i in range(0, n_chains):

        chain_atom23_pos = atom23_pos[chain_indices == i]
        #chain_atom3_pos = chain_atom23_pos[:, [0, 5, 12], :]

        chain_res_types = [res_types[index] for index in range(len(res_types)) if chain_indices[index] == i]

        chain_atom3_pos = atom3_pos[chain_indices == i]

        # Use raw structure
        #chain_atom3_pos = raw_atom23_pos[raw_chain_indices == i][:, [0, 5, 12], :]

        r_atom23_pos, r_atom23_mask = full_atom(
            chain_atom3_pos, 
            chain_res_types, 
            lib_dir=lib_dir, flag="", 
            natype='dna' if final_chains_is_dna[i] else 'rna',
        )

        chains_atom23_pos.append(r_atom23_pos)
        chains_atom23_mask.append(r_atom23_mask)
        chains_is_dna.append(final_chains_is_dna[i])
        chains_res_types.append(chain_res_types)


    # TODO add confidence prediction as bfactor
    pass

    # Write final
    chains_bfactors = None
    #fpdbout = pjoin(fout, "reassign_helix.pdb")
    fpdbout = pjoin(fout, "reassign_helix.cif")
    chains_atom23_to_pdb(
        fpdbout,
        chains_atom23_pos=chains_atom23_pos,
        chains_atom23_mask=chains_atom23_mask,
        chains_res_types=chains_res_types,
        chains_bfactors=chains_bfactors,
        chains_is_dna=chains_is_dna,
    )
    print("# Output final chains to {}".format(fpdbout))

    tend = time.time()
    print("# Time consuming {:.4f}".format(tend - tstart))

