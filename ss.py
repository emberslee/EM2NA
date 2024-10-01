import os
import sys
import math
import tempfile
import argparse
import subprocess
import numpy as np
np.random.seed(42)
from pdbio import chains_atom23_to_pdb, nuc_type_index_to_1, read_atom23_na, nuc_type_index_to_1, nuc_type_1_to_index
from opt_seq_assign_helix import get_ss_by_cssr
from seqio import read_secstr, read_fasta, ss_to_pairs
from ortools.sat.python import cp_model
from seqio import align
from opt_seq_assign import seq_identity
from grid import get_grid_value
from utils import parse_map


patterns = [
    ('.', '.'), 
    ('(', ')'), ('[', ']'), ('{', '}'), ('<', '>'),
    ('A', 'a'), ('B', 'b'), ('C', 'c'), ('D', 'd'),
    ('E', 'e'), ('F', 'f'), ('G', 'h'), ('H', 'h'),
    ('I', 'i'), ('J', 'j'), ('K', 'k'), ('L', 'l'),
    ('M', 'm'), ('N', 'n'), ('O', 'o'), ('P', 'p'),
    ('Q', 'q'), ('R', 'r'), ('S', 's'), ('T', 't'),
    ('U', 'u'), ('V', 'v'), ('W', 'w'), ('X', 'x'),
    ('Y', 'y'), ('Z', 'z'),
]
patterns_left = []
patterns_right = []
patterns_left_to_idx = {}
patterns_right_to_idx = {}
patterns_idx_to_left = {}
patterns_idx_to_right = {}
for i, p in enumerate(patterns):
    l = p[0]
    r = p[1]
    patterns_left.append(l)
    patterns_right.append(r)
    patterns_left_to_idx[l] = i
    patterns_right_to_idx[r] = i
    patterns_idx_to_left[i] = l
    patterns_idx_to_right[i] = r

patterns_left_to_right = {}
for p in patterns:
    l = p[0]
    r = p[1]
    patterns_left_to_right[l] = r

patterns_right_to_left = {}
for l, r in patterns_left_to_right.items():
    patterns_right_to_left[r] = l


def get_level(ss, l):
    if ss[l[0]] != '.' or ss[l[1]] != '.':
        return 0
    for level in range(1, len(patterns)):
        score = 0
        flag = 1
        for i in range(l[0]+1, l[1]):
            if ss[i] == patterns[level][0]:
                score += 1
            elif ss[i] == patterns[level][1]:
                score -= 1
            if score < 0:
                flag = 0
                break
        if score != 0:
            flag = 0
        if flag == 1:
            return level
    return 0

def pairs_to_ss(l, pairs):
    size = l
    ss = ['.'] * l
    for l in pairs:
        level = get_level(ss, l)
        if level != 0:
            ss[l[0]] = patterns[level][0]
            ss[l[1]] = patterns[level][1]
    return "".join(ss)

def filter_pairs(pairs):
    ret = []
    idx_max = 0
    for p in pairs:
        idx_max = max(idx_max, p[0])
        idx_max = max(idx_max, p[1])
    length = idx_max + 1
    flag = [0] * length
    for p in pairs:
        flag[p[0]] = 1
        flag[p[1]] = 1

    for p in pairs:
        l, r = p
        if l - 1 >= 0 and flag[l - 1] == 0 and \
            l + 1 < length and flag[l + 1] == 0:
            continue
        if r - 1 >= 0 and flag[r - 1] == 0 and \
            r + 1 < length and flag[r + 1] == 0:
            continue
        ret.append(p)
    return ret

def ss_to_sse(ss, n_pair_min=2, verbose=False):
    sses = []
    l = len(ss)
    stacks = [[] for i in range(len(patterns))]
    i = 0
    while i < l:
        ch = ss[i]
        if ch in patterns_left_to_idx.keys():
            idx = patterns_left_to_idx[ch]
            stacks[idx].append((i, ch))
        # is a right bracket
        else:
            k = i
            while k < l and ss[k] == ch:
                k += 1
            n = k - i
            print("# Intend to push out {} pairs".format(n))
            # push out n chs
            idx = patterns_right_to_idx[ch]
            head = []
            while n > 0:
                if not stacks[idx]:
                    raise Exception("Error invalid secondary structure -> {}".format(ss))
                top = stacks[idx][-1]
                head.append(top)
                stacks[idx].pop()
                n -= 1
            head = head[::-1]

            tail = []
            for m in range(i, k):
                tail.append((m, ss[m]))

            sse = head + tail
            sses.append(sse)
            i = k - 1
        i += 1

    # the left stacks should be empty
    valid = True
    for i, s in enumerate(stacks):
        if i == 0:
            continue
        if verbose:
            print("# Check stack for idx {} -> {}".format(i, s))
        if s:
            raise Exception("Error there are unmatched brackets in the stack maybe you have input invalid secondary structure -> {}".format(ss))
    sses = [sse for sse in sses if len(sse) >= 2 * n_pair_min]
    print("# Split to {} SSEs".format(len(sses)))
    return sses 

def print_sse(sse):
    idx_str = "# "
    bracket_str = "# "
    for site in sse:
        idx_str += "{} ".format(site[0])
        bracket_str += "{}".format(site[1])
    print(idx_str)
    print(bracket_str)


def main(args):
    fin = args.i
    lib_dir = args.l
    if lib_dir is None:
        lib_dir = os.path.abspath(os.path.dirname(__file__))
    lib_dir = os.path.abspath(os.path.expanduser(lib_dir))

    # parse input seq
    fseq = args.seq
    seqs = read_fasta(fseq)
    if len(seqs) == 1:
        seq = seqs[0]
    elif len(seqs) >= 2:
        raise Exception("Error found >=2 sequences as input can only handle 1 as input")
    else:
        raise Exception("Error cannot find any sequence in {}".format(fseq))
    tgt_seq = seq

    # parse input sec
    fsec = args.sec
    secs = read_secstr(fsec)
    if len(secs) == 1:
        sec = secs[0]
    elif len(secs) >= 2:
        raise Exception("Error found >=2 secondary structures as input can only handle 1 as input")
    else:
        raise Exception("Error cannot find any secondary structure in {}".format(fsec))

    # convert to pair repr.
    tgt_pairs = ss_to_pairs(sec)
    tgt_pairs = filter_pairs(tgt_pairs)
    tgt_sec = pairs_to_ss(l=len(sec), pairs=tgt_pairs)
    print(tgt_pairs)

    print("# Target sequence")
    print(seq)
    print("# Target secondary structure (initial)")
    print(sec)
    print("# Target secondary structure (transformed)")
    print(tgt_sec)

    # split to sses
    n_pair_min = 2
    tgt_sses = ss_to_sse(tgt_sec, n_pair_min=n_pair_min)
    for sse in tgt_sses:
        print_sse(sse)

    if not len(seq) == len(sec):
        raise Exception("# Error length of sequence and secondary structure not equal")


    # read input structure
    atom_pos, atom_mask, res_type, _, _ = read_atom23_na(fin)
    tgt_initial_full_seq = [nuc_type_index_to_1[x] for x in res_type]

    # parse
    temp_dir = './temp'
    pairs_dir = os.path.join(temp_dir, "pairs.txt")
    if not os.path.exists(pairs_dir):
        print("# Parsing sec of build model")
        pairs = get_ss_by_cssr(atom_pos, lib_dir=lib_dir, out_dir=temp_dir)
        pairs = filter_pairs(pairs)
        print(pairs)
        with open(pairs_dir, "w") as f:
            for p in pairs:
                f.write("{} {}\n".format(p[0], p[1]))
    else:
        print("# Read cache")
        pairs = []
        with open(pairs_dir, "r") as f:
            for line in f.readlines():
                p0, p1 = line.strip().split()
                pairs.append((int(p0), int(p1)))

    ss = pairs_to_ss(l=len(atom_pos), pairs=pairs)
    print("# Query secondary structure")
    print(ss)
    cov = 2 * len(pairs) / (len(atom_pos) + 1e-6)
    print("# Detected ss covers {:.4f} % residues".format(cov * 100.))

    # 1. split build model into substructure
    # split to sses
    sses = ss_to_sse(ss, n_pair_min=n_pair_min)
    for sse in sses:
        print_sse(sse)

    # 2. align substructure using IP
    # 2.1 calculate alignment score for every two substrctures
    # target: ss from build structure
    # query: user input ss
    m = len(tgt_sses)
    n = len(sses)

    # TODO
    # setup up scoring matrix
    substr_score_mtx = np.zeros((m, n), dtype=np.float32)
    rev_score_mtx = np.zeros((m, n, 2), dtype=np.float32)
    all_alignment = [[[None, None] for k in range(n)] for i in range(m)]
    for i in range(m):
        # get index of build structure
        tgt_sse = tgt_sses[i]
        tgt_idxs = [x[0] for x in tgt_sse]
        #print(tgt_idxs)

        # get res type or logits
        tgt_seq = [nuc_type_index_to_1[x] for x in res_type]
        tgt_res_type = [tgt_seq[x] for x in tgt_idxs]
        tgt_res_type = "".join(tgt_res_type)

        for k in range(n):
            # get index of query sequence and ss
            qry_sse = sses[k]
            qry_idxs_left = [x[0] for x in qry_sse if x[1] in patterns_left]
            qry_idxs_right = [x[0] for x in qry_sse if x[1] in patterns_right]
            qry_idxs = qry_idxs_left + qry_idxs_right

            # get res type
            qry_res_type = [seq[x] for x in qry_idxs]
            qry_res_type = "".join(qry_res_type)

            # reverse restype
            qry_idxs_rev = qry_idxs_right + qry_idxs_left
            qry_res_type_rev = [seq[x] for x in qry_idxs_rev]
            qry_res_type_rev = "".join(qry_res_type_rev)

            #print(tgt_res_type)
            #print(qry_res_type)
            #print(qry_res_type_rev)

            #alignment     = align(tgt_res_type, qry_res_type,     gap_open=-0.5, gap_extend=-0.5)[0]
            #alignment_rev = align(tgt_res_type, qry_res_type_rev, gap_open=-0.5, gap_extend=-0.5)[0]

            # 0722
            # align only half
            alignment     = align(tgt_res_type, qry_res_type)[0]
            alignment_rev = align(tgt_res_type, qry_res_type_rev)[0]
            #print("# len {} {}".format(len(tgt_res_type), len(qry_res_type)))
            #print("# forward {:.4f}".format(alignment.score))
            #print("# bckward {:.4f}".format(alignment_rev.score))

            rev_score_mtx[i][k][0] = alignment.score
            rev_score_mtx[i][k][1] = alignment_rev.score

            all_alignment[i][k][0] = alignment
            all_alignment[i][k][1] = alignment_rev

            substr_score_mtx[i][k] = max(alignment.score, alignment_rev.score)
    #exit()

    # when matching tgt_sec to sec
    # tgt_sec should be also reversed to match
    # e.g.
    # if the build structure is presented (in .cif) as:
    # AGCUGCU ->
    # |||||||  |
    # UCGACGA <-
    # and if the target seq and sec are presented as:
    # AGCUGCU -> (part1)
    #  ||||||  |
    # ACGACGA <- (part2)
    # they are nearly perfectly matched
    # AGCUGCUAGCAGCU
    # :::::::::::::.
    # AGCUGCUAGCAGCA
    # and if the target seq and sec are presentaed as:
    # AGCAGCA -> (reversed part2)
    #  ||||||  |
    # UCGUCGA <- (reversed part1)
    # the match becomes
    # AGCUGCUAGCAGCU
    # :::.::.:::.:: 
    # AGCAGCAAGCUGCA
    # which is worse than initial match
    # but they SHOULD be matched logically
    # so we should either transform build seq/sec or input seq/sec to align twice
    

    # rescale and covert to integer
    rescale_substr_score_mtx = (substr_score_mtx * 10).astype(np.int32)

    # add model
    model = cp_model.CpModel()

    # add decision variable - correspondence
    x = [[model.NewBoolVar(f'x[{i}][{k}]') for k in range(n)] for i in range(m)]

    # TODO add sequence restraint
    # 1. abs(len(i), len(m)) < 2 * max_pair
    max_pair = 3
    for i in range(m):
        for k in range(n):
            if abs(len(tgt_sses[i]) - len(sses[k])) > 2 * max_pair:
                model.Add(x[i][k] == 0)

    # add row/col restraint
    for i in range(m):
        model.Add(sum(x[i][k] for k in range(n)) <= 1)

    for k in range(n):
        model.Add(sum(x[i][k] for i in range(m)) <= 1)

    # set objective
    objectives = []
    for i in range(m):
        for k in range(n):
            objectives.append(
                x[i][k] * rescale_substr_score_mtx[i][k]
            )
    model.Maximize(sum(objectives))
 
    # solve
    time_limit = 300
    num_workers = 2
    log_search_progress = True

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = num_workers
    solver.parameters.log_search_progress = log_search_progress
    solution_printer = cp_model.ObjectiveSolutionPrinter()
    status = solver.SolveWithSolutionCallback(model, solution_printer)

    corr = np.zeros((m, n), dtype=np.int32)
    if (status == cp_model.OPTIMAL) or \
        (status == cp_model.FEASIBLE) or \
        (status == cp_model.INTERRUPTED and solver.ResponseStats().has_feasible_solution()):
        for i in range(m):
            for k in range(n):
                corr[i][k] = 1 if solver.BooleanValue(x[i][k]) == 1 else 0
    #print(corr)
    print(corr.sum())

    # read input aamap
    faamap = args.aamap
    aalogits0, aaorigin0, _, aavsize0 = parse_map(faamap, False, None)
    print("# Read map from {}".format(faamap))

    # Recalculate aaprobs
    faaprob = args.aaprob
    aaprob = np.load(faaprob)
    aaprob = aaprob['prob']
    print("# done reading aaprob from {}".format(faaprob))

    aaprobs = []
    npos = atom_pos[:, 12]
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
    print("# done assign aaprob for each nucleotide")

    # 3. postprocess
    # when matched sequence have different length
    # if shorter, extend both side
    #keep = [False] * len(atom_pos)
    keep = [True] * len(atom_pos)
    new_res_type = res_type.copy()
    for i in range(m):
        k = np.argmax(corr[i])
        #print('-'*50)
        #print( len(tgt_sses[i]) )
        #print( len(sses[k]) )
        #print(substr_score_mtx[i][k])
        #print(rev_score_mtx[i][k][0], rev_score_mtx[i][k][1])
        if corr[i][k] != 1:
            continue

        if np.abs(substr_score_mtx[i][k] - rev_score_mtx[i][k][0]) < 1e-3:
            sel = 0
        else:
            sel = 1
        #print(sel)

        # get alignment
        alignment = all_alignment[i][k][sel]
        #print(alignment)

        seqA, seqB = alignment[:2]
        aligned_seq = []
        for l in range(len(seqA)):
            if seqA[l] != '-':
                aligned_seq.append(seqB[l] if seqB[l] != '-' else seqA[l])
       
        tgt_aligned_idxs = [x[0] for x in tgt_sses[i]]
        tgt_aligned_seq = "".join(aligned_seq)
        
        #print(len(tgt_aligned_idxs), len(tgt_aligned_seq))
        assert len(tgt_aligned_idxs) == len(tgt_aligned_seq)
        print("# before corretiong bp seq = {}".format(tgt_aligned_seq))
        #l = len(tgt_aligned_seq) // 2
        #print(tgt_aligned_seq[:l])
        #print(tgt_aligned_seq[l:][::-1])
      
        tgt_initial_seq = [nuc_type_index_to_1[res_type[x]] for x in tgt_aligned_idxs]
        print(tgt_initial_seq)
        tgt_initial_seq = "".join(tgt_initial_seq)
        print("# initial sequence is      = {}".format(tgt_initial_seq))

        # correct_bp
        if args.bp_correct:
            bp_res_type = [x for x in tgt_aligned_seq]
            watson_crick_pair = {
                'A': 'U',
                'U': 'A',
                'C': 'G',
                'G': 'C',
            }
            n_aligned = len(tgt_aligned_idxs)

            probs_aligned = []
            probs_initial = []
            for l in range(n_aligned // 2):
                base0 = tgt_aligned_seq[l]
                base1 = tgt_aligned_seq[n_aligned - l - 1]
                base0 = nuc_type_1_to_index[base0]
                base1 = nuc_type_1_to_index[base1]

                # 1. use assigned base type
                idx0 = tgt_aligned_idxs[l]
                idx1 = tgt_aligned_idxs[n_aligned - l - 1]

                prob0 = aaprobs[idx0][base0]
                prob1 = aaprobs[idx1][base1]
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
    
                bp_res_type[l]                 = rbase0
                bp_res_type[n_aligned - l - 1] = rbase1

                idx0 = tgt_aligned_idxs[l]
                # calculate probs of matched
                probs_aligned.append(aaprobs[idx0][nuc_type_1_to_index[ tgt_aligned_seq[l] ]])
                # calculate probs of initial
                probs_initial.append(aaprobs[idx0][nuc_type_1_to_index[ tgt_initial_seq[l] ]])

            probs_aligned = np.mean(probs_aligned)
            probs_initial = np.mean(probs_initial)

            print(True if probs_aligned > probs_initial else False)
            print(probs_aligned, probs_initial)

            if probs_aligned > probs_initial:
                tgt_aligned_seq = [nuc_type_index_to_1[x] for x in bp_res_type]
                tgt_aligned_seq = "".join(tgt_aligned_seq)
            else:
                tgt_aligned_seq = tgt_initial_seq
            print("# after  corretiong bp seq = {}".format(tgt_aligned_seq))
 

        # for each aligned site
        for l in range(len(tgt_aligned_idxs)):
            idx = tgt_aligned_idxs[l]
            new_res_type[idx] = nuc_type_1_to_index[tgt_aligned_seq[l]]
            keep[idx] = True

    new_atom_mask = np.zeros_like(atom_mask, dtype=bool)
    new_atom_mask[..., [0, 5, 12]] = True

    chains_atom23_to_pdb(
        "test.cif",
        chains_atom23_pos=[atom_pos[keep, ...]],
        chains_atom23_mask=[new_atom_mask[keep, ...]],
        chains_res_types=[[nuc_type_index_to_1[x] for x in new_res_type[keep, ...]]],
    )

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="Input RNA structure")
    parser.add_argument("-l", help="Lib dir")
    parser.add_argument("--seq", "-seq", help="Input sequence")
    parser.add_argument("--sec", "-sec", help="Input secondary structure")
    parser.add_argument("--tgt_sec", help="Target secondary structure")
    parser.add_argument("--aamap", "-aamap", help="Input aa map", required=True)
    parser.add_argument("--aaprob", "-aaprob", help="Inpu aaprob", required=True)
    parser.add_argument("--bp_correct", action='store_true')
    #
    parser.add_argument("-t", help="Time limit", type=int)
    parser.add_argument("-v", help="Verbose", action='store_true')
    parser.add_argument("-nt", help="Num threads", type=int)
    args = parser.parse_args()
    main(args)
