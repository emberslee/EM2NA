from Bio import pairwise2
from Bio.Seq import Seq
import warnings
warnings.filterwarnings("ignore")
import numpy as np


def readlines(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    lines = [line.strip() for line in lines if line.strip()]
    return lines

def read_lines(filename):
    lines = readlines(filename)
    return lines

def read_fasta(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    seqs = []
    seq = ""
    for line in lines:
        if line.startswith('>'):
            if len(seq) > 0:
                seqs.append(seq)
            seq = ""
            continue
        seq += line.strip()
    if len(seq) > 0:
        seqs.append(seq)
    return seqs

def write_lines_to_file(lines, filename):
    with open(filename, 'w') as f:
        for line in lines:
            f.write(line.strip())
            f.write("\n")

def write_seq_to_file(seq, filename):
    lines = [">seq"]
    lines.append(seq.strip())
    write_lines_to_file(lines, filename)

def read_secstr(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    seqs = []
    seq = ""
    for line in lines:
        if line.startswith('>'):
            if len(seq) > 0:
                seqs.append(seq)
            seq = ""
            continue
        seq += line.strip()
    if len(seq) > 0:
        seqs.append(seq)
    return seqs
    '''
    with open(filename, 'r') as f:
        lines = f.readlines()
    secs = []
    for line in lines:
        if len(line.strip()) > 0:
            secs.append(line.strip())
    return secs
    '''

def toupper(seq):
    return seq.upper()

def tolower(seq):
    return seq.lower()

def is_na(seq):
    nucs = ['A', 'G', 'C', 'U', 'T', 'I', 'N', 'X', 'a', 'g', 'c', 'u', 't', 'i', 'n', 'x']
    for s in seq:
        if s not in nucs:
            return False
    return True

def is_dna(seq):
    return is_na(seq) and ('T' in seq or 't' in seq)

def is_rna(seq):
    return is_na(seq) and 'T' not in seq and 't' not in seq

def align(a, b, mode='global', match=2, mismatch=-1, gap_open=-0.5, gap_extend=-0.1, n_classes=4, setting=1, **kwargs):
    seq1 = Seq(a)
    seq2 = Seq(b)
    # Get match dict
    #print("# Use setting = {}".format(setting))
    if setting == 0:
        match_dict = {
            ('A', 'A'): 2, ('A', 'G'): 1, ('A', 'C'):-1, ('A', 'U'):-1, ('A', 'T'):-1, 
            ('G', 'A'): 1, ('G', 'G'): 2, ('G', 'C'):-1, ('G', 'U'):-1, ('G', 'T'):-1, 
            ('C', 'A'):-1, ('C', 'G'):-1, ('C', 'C'): 2, ('C', 'U'): 1, ('C', 'T'): 1, 
            ('U', 'A'):-1, ('U', 'G'):-1, ('U', 'C'): 1, ('U', 'U'): 2, ('U', 'T'): 1, 
            ('T', 'A'):-1, ('T', 'G'):-1, ('T', 'C'): 1, ('T', 'U'): 1, ('T', 'T'): 2, 
        }
    else:
        """
        aa_match_matrix = (np.asarray(
            [
                [0.5312, 0.2079, 0.1467, 0.1141], 
                [0.1871, 0.5935, 0.1475, 0.0719], 
                [0.1258, 0.1565, 0.5934, 0.1243], 
                [0.1673, 0.1445, 0.2536, 0.4346], 
            ],
            dtype=np.float32
        ) - 0.10) * 2
        """


        """
        aa_match_matrix = (np.array(
            [
                [0.5000, 0.3000, 0.1000, 0.1000],
                [0.2500, 0.6000, 0.1000, 0.0500],
                [0.1500, 0.1500, 0.6000, 0.2000],
                [0.1500, 0.1000, 0.3000, 0.4500],
            ]
        ) - 0.10) * 4
        """


        # 2024-01-24
        """
        465 310 167 194
        303 569 177 147
        160 161 477 252
        194 121 365 407
        """
        aa_match_matrix = np.array(
            [[0.20466549, 0.13293310, 0.07639524, 0.08726946],
             [0.12993139, 0.23787625, 0.07880677, 0.06438896],
             [0.07319305, 0.07168299, 0.22714286, 0.11792232],
             [0.08726946, 0.05300044, 0.17080019, 0.18721251]]
        )
        aa_match_matrix = (aa_match_matrix - 0.05) * 10

        """
        aa_match_matrix = np.asarray(
            [[0.21395669, 0.07300743, 0.02529714, 0.03512650],
             [0.07055529, 0.28289710, 0.03002405, 0.02208786],
             [0.02288233, 0.02619477, 0.23350690, 0.05990050],
             [0.03239937, 0.01781308, 0.09237735, 0.19725757]]
        )
        aa_match_matrix = (aa_match_matrix - 0.05) * 10
        """

        names = ['A', 'G', 'C', 'U', 'T']
        match_dict = dict()
        for i in range(4):
            for k in range(4):
                match_dict[(names[i], names[k])] = aa_match_matrix[i][k]


    if mode == 'global':
        alignments = pairwise2.align.globalds(seq1, seq2, match_dict, gap_open, gap_extend, **kwargs)
    elif mode == 'local':
        alignments =  pairwise2.align.localds(seq1, seq2, match_dict, gap_open, gap_extend, **kwargs)
    else:
        raise ValueError
    return alignments


def multi_align(a, b, mode='global', match=2, mismatch=-1, gap_open=-0.5, gap_extend=-0.1, n_classes=4, top=20):
    seq1 = Seq(a)
    seq2 = Seq(b)
    # Get match dict
    match_dict = {
        ('A', 'A'): 2, ('A', 'G'): 1, ('A', 'C'):-1, ('A', 'U'):-1, ('A', 'T'):-1, 
        ('G', 'A'): 1, ('G', 'G'): 2, ('G', 'C'):-1, ('G', 'U'):-1, ('G', 'T'):-1, 
        ('C', 'A'):-1, ('C', 'G'):-1, ('C', 'C'): 2, ('C', 'U'): 1, ('C', 'T'): 1, 
        ('U', 'A'):-1, ('U', 'G'):-1, ('U', 'C'): 1, ('U', 'U'): 2, ('U', 'T'): 1, 
        ('T', 'A'):-1, ('T', 'G'):-1, ('T', 'C'): 1, ('T', 'U'): 1, ('T', 'T'): 2, 
    }
    if mode == 'global':
        alignments = pairwise2.align.globalds(seq1, seq2, match_dict, gap_open, gap_extend, one_alignment_only=False)
    elif mode == 'local':
        alignments =  pairwise2.align.localds(seq1, seq2, match_dict, gap_open, gap_extend, one_alignment_only=False)
    else:
        raise ValueError
    print(len(alignments))
    return alignments[:top] if top > 1 else [alignments[top]]



def read_era_result(filename):
    result = readlines(filename)
    if len(result) < 5:
        raise "ERA result is not valid"
    return result[-5:]


def read_usalign_result(filename):
    result = readlines(filename)
    if len(result) < 3:
        raise "USalign result is not valid"
    return result[-4:-1]



'''
    input
    seqA : ---AAAGASG-ASGASG--ASFASF--
    seqB : ASDASO--S-ASF--ASFAFSF-ASAD
    return
    seqA : AAAGASG-ASGASG--ASFASF
    seqB : ASO--S-ASF--ASFAFSF-AS
'''
def remove_bar(seq1, seq2, sec2=None):
    assert len(seq1) == len(seq2)
    if sec2 is not None:
        assert len(seq2) == len(sec2)

    # Remove bar both terminus in seq1
    n = len(seq1)
    ret_seq1 = []
    ret_seq2 = []
    ret_sec2 = []
    for i in range(n):
        if seq1[i] != '-':
            ret_seq1.append(seq1[i])
            ret_seq2.append(seq2[i])
            if sec2 is not None:
                ret_sec2.append(sec2[i])

    # s1 does not have '-'
    # s2 may have '-'
    n = len(ret_seq2)
    for i in range(n):
        if ret_seq2[i] == '-':
            ret_seq2[i] = ret_seq1[i]
            if sec2 is not None:
                ret_sec2[i] = '.'

    ret_seq1 = "".join(ret_seq1)
    ret_seq2 = "".join(ret_seq2)
    ret_sec2 = "".join(ret_sec2)
    if sec2 is not None:
        return ret_seq1, ret_seq2, ret_sec2
    else:
        return ret_seq1, ret_seq2


# Simple removes '-' in seq2
def remove_bar_simple(seq1, seq2, sec2=None):
    assert len(seq1) == len(seq2)
    if sec2 is not None:
        assert len(seq2) == len(sec2)

    ret_seq1 = "".join([x for x in seq1 if x != '-'])

    ret_seq2 = []
    ret_sec2 = []

    for i in range(len(seq2)):
        if seq2[i] == '-':
            continue
        ret_seq2.append(seq2[i])
        if sec2 is not None:
            ret_sec2.append(sec2[i])

    ret_seq2 = "".join(ret_seq2)
    ret_sec2 = "".join(ret_sec2)

    if sec2 is not None:
        return ret_seq1, ret_seq2, ret_sec2
    else:
        return ret_seq1, ret_seq2


# If seq1[i] != '-' and seq2[i] != '-'
# seq1[i] = seq2[i]
def remove_bar_match(seq1, seq2):
    assert len(seq1) == len(seq1)
    seqA = []
    seqB = []
    for i in range(len(seq1)):
        if seq1[i] != '-' and seq2[i] != '-':
            seqA.append(seq1[i])
            seqB.append(seq2[i])
        else:
            seqA.append(seq1[i])
            seqB.append(seq1[i])

    seqA = "".join(seqA)
    seqB = "".join(seqB)

    seqA = seqA.replace('-', '')
    seqB = seqB.replace('-', '')

    return seqA, seqB


# Parse a second structure
def ss_to_pairs(ss):
    ret = []
    s = []

    for i in range(len(ss)):
        if ss[i] == '.':
            continue

        # Find bracket
        if ss[i] in ['{', '[', '<', '(']:
            s.append( (ss[i], i) )
            continue
        else:
            symbol = '.';
            if ss[i] == ')':
                symbol = '('
            elif ss[i] == ']':
                symbol = '['
            elif ss[i] == '}':
                symbol = '{'
            elif ss[i] == '>':
                symbol = '<'

            # Find nearest left bracket
            s0 = []
            while s:
                ch = s[-1][0]
                ii = s[-1][1]
                if ch == symbol:
                    ret.append( (ii, i) )
                    s.pop()
                    break
                else:
                    s0.append( s[-1])
                    s.pop()

            # Recover s
            while s0:
                s.append( s0[-1] )
                s0.pop()

    # Sort by first number
    ret.sort(key=lambda x:x[0])
    return ret


def align_seq_and_sec(seq, sec):
    assert len(seq) >= len(sec)
    ret = []
    i, k = 0, 0
    while i < len(seq):
        if seq[i] != '-':
            ret.append(sec[k])
            i+=1
            k+=1
        else:
            ret.append('&')
            i+=1
    # Postprocess
    if len(ret) > len(seq):
        ret = ret[:len(seq)]
    elif len(ret) < len(seq):
        ret.extend(['&']*(len(seq)-len(ret)))
    return "".join(ret)


def convert_to_rna_seq(seq):
    return seq.replace('T', 'U')

def convert_to_dna_seq(seq):
    return seq.replace('U', 'T')


def get_multi_seq_align_among_candidates(seq1, seqs, clean_gap=True):
    # Determine which sequence best fits the segment
    score0 = -1e6
    seqA0 = None
    seqB0 = None
    k0 = None

    for k, seq2 in enumerate(seqs):
        # Target sequence is >= query sequence
        '''
        if len(seq2) < len(seq1):
            continue
        '''

        alignments = multi_align(seq1, seq2)
        for i, alignment in enumerate(alignments):
            seqA, seqB, score, start, end = alignment[:5]
            print(seqA)
            print(seqB)
            print(score)

    seqA = seqA0
    seqB = seqB0






def get_seq_align_among_candidates(seq1, seqs, **kwargs):
    # Determine which sequence best fits the segment
    score0 = -1e6
    seqA0 = None
    seqB0 = None
    k0 = None

    for k, seq2 in enumerate(seqs):
        # Target sequence is >= query sequence
        '''
        if len(seq2) < len(seq1):
            continue
        '''

        alignment = align(seq1, seq2, gap_open=-1.0, gap_extend=-0.5, **kwargs)[0]
        seqA, seqB, score, start, end = alignment[:5]

        if score > score0:
            score0 = score
            seqA0 = seqA[start:end+1]
            seqB0 = seqB[start:end+1]
            k0 = k

    seqA = seqA0
    seqB = seqB0
    #print(seqA)
    #print(seqB)

    # In case that given sequence is invalid
    if seqA is None or \
       seqB is None:
        seqA = seq1
        seqB = seq1

    seqA, seqB = remove_bar(seqA, seqB)
    #seqA, seqB = remove_bar_match(seqA, seqB)

    # In case that seqB is shorted
    seqB = seqB + 'U'*(len(seq1)-len(seqB))
    return seqA, seqB, score0, k0


def format_seqs(seqs):
    # 0 Convert to upper
    seqs = [seq.upper() for seq in seqs]

    # 1. Convert chars not in "AGCUT" to "AGCUT"
    seqs0 = []
    for seq in seqs:
        seq0 = ""
        for s in seq:
            if s in ['A', 'G', 'C', 'U', 'T']:
                seq0 += s
            else:
                seq0 += 'U'
        seqs0.append(seq0)
    return seqs0


# some utils for ranges
def is_valid_interval(length, interval):
    start, end = interval
    return 0 <= start < length and 0 <= end < length and start < end

def valid_interval(length, interval):
    start, end = interval
    if start < 0:
        start = 0
    if end > length:
        end = length
    return (start, end)

def remove_intervals(length, intervals):
    intervals = [valid_interval(length, x) for x in intervals]
    flag = [True] * length
    for interval in intervals:
        start, end = interval
        for k in range(start, end):
            flag[k] = False
    ret = []
    k = 0
    while k < length:
        start = k
        interval = []
        while start < length and not flag[start]:
            start += 1
        while start < length and flag[start]:
            interval.append(start)
            start += 1
        if len(interval) >= 1:
            ret.append([interval[0], interval[-1]])
        k = start
    return ret


if __name__ == '__main__':
    ss = "(((((((((((....))))...)))))))"
    print(ss)
    print(ss_to_pairs(ss))

    ss = "...((((.....)))))))))...)))))"
    print(ss)
    print(ss_to_pairs(ss))

    ss =  ".................."
    print(ss)
    print(ss_to_pairs(ss))

    ss = "<<<<<...((((((..((((.....>>>>>>>....)))))))))))))))))))))))))))"
    print(ss)
    print(ss_to_pairs(ss))

    sec = align_seq_and_sec("ASFASF---SFA-ASF--", "(((...))))...")
    print("ASFASF---SFA-ASF--")
    print(sec)



    s1 = "GGUCCCC---GCCCUU-CGGCCGCC-CC-CCGGCA----AA-AGAACCGCCUGGAAAAAAAAAGACGGG-GGGA-UUAGCCC-CGUGGCUCCG-GCCACCGUAAGGGGUUGG--GAGGCCCCUUUU--UA--CUACCGCCUU---AUGGCGAAAAGG--GGCCG-GGAAGUUU-UAAG-CCUAACU-CCGGAAAACC-GCCU-UGGGGCCCCGUUUUA----CUC-CCCUC---GGGUUUAAGGCGGUU-U-"
    s2 = "GGU---UAAAG-CCUUAUGGUCGCUACCAUU-GCACUCCGGUAG---CG-UU--AAAAGGGAAGACGGGUGAGAAU---CCCGCGCAGCCCC-CGCUACUGUGAGGGA--GGACGAAGCCC----UAGUAAGCCACUGCC--GAAA-GGUGGGAAGGCAGG--GUGGAGG---AUGAGUCCCGA--GCCAGGAGACCUGCC-AUAAG-----GUUUUAGAAGUUCGCCUUCGGGGGG---AAGGUGA--ACA"

    s1, s2 = remove_bar_simple(s1, s2)
    print(s1)
    print(s2)


    
    ranges = remove_intervals(128, [[4, 5], [4, 7], [8, 15]])
    print(ranges)
