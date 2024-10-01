import os
import numpy as np
from pdbio import read_atom23_na

def pairwise_distances(a, b):
    """
    a, b : np.ndarray of shape [L, d], [M, d]
    """
    assert a.ndim == b.ndim
    assert a.shape[-1] == b.shape[-1]
    distances = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1)
    return distances


def kabsch(P, Q):
    P = np.asarray(P, dtype=np.float32)
    Q = np.asarray(Q, dtype=np.float32)
    assert len(P) == len(Q)
    assert len(P) > 0

    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)
    P_centered = P - centroid_P
    Q_centered = Q - centroid_Q
    H = P_centered.T.dot(Q_centered)
    U, S, VT = np.linalg.svd(H)
    R = U.dot(VT).T
    if np.linalg.det(R) < 0:
        VT[2,:] *= -1
        R = U.dot(VT).T
    t = centroid_Q - R.dot(centroid_P)

    return  R, t

def kabsch_rmsd(P, Q):
    P = np.asarray(P, dtype=np.float32)
    Q = np.asarray(Q, dtype=np.float32)

    R, T = kabsch(P, Q)
    rP = P @ R.T + T
    
    return np.sqrt(np.sum((rP - Q) ** 2))



def helix_form(atom3_pos, pairs, lib_dir=None):
    # At least two pairs
    #if len(pairs) < 3:
    #    return 0

    assert lib_dir is not None, "# Error Please specify lib directory"
    fdna_template = os.path.join(lib_dir, "lib", "DNA", "ACGT.dna.pdb")
    frna_template = os.path.join(lib_dir, "lib", "RNA", "ACGU.rna.pdb")
    dna_atom23_pos, _, _, _, _, _ = read_atom23_na(fdna_template, check_dna=True)
    rna_atom23_pos, _, _, _, _, _ = read_atom23_na(frna_template, check_dna=True)

    dna_atom3_pos = dna_atom23_pos[:, [0, 5, 12], :]
    rna_atom3_pos = rna_atom23_pos[:, [0, 5, 12], :]

    dna_crds = dna_atom3_pos[[0, 7], :].reshape(-1, 3) # (6, 3)
    rna_crds = rna_atom3_pos[[0, 7], :].reshape(-1, 3) # (6, 3)

    #print(len(dna_crds))
    #print(len(rna_crds))

    forms = []
    # 0. Using template base-pair
    for p in pairs:
        query_pos = atom3_pos[[p[0], p[1]]].reshape(-1, 3)

        #print(dna_crds.shape)
        #print(query_pos.shape)

        dna_rmsd = kabsch_rmsd(dna_crds, query_pos)
        rna_rmsd = kabsch_rmsd(rna_crds, query_pos)
  
        #print("{:.4f} {:.4f}".format(dna_rmsd, rna_rmsd))
 
        forms.append(0 if rna_rmsd < dna_rmsd else 1)

    # 1. Using geometry
    # todo

    return forms


if __name__ == '__main__':
    a = np.array([[1,4,7], [2,5,8]])
    b = np.array([[1,2,1], [2,5,2]])

    R, T = kabsch(a, b)
    print(R)
    print(T)
    print(kabsch_rmsd(a, b))


    c = a @ R.T + T
    
    print(b)
    print(c)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    print(script_dir)   

    query_atom23_pos, _, _, _, _, _ = read_atom23_na(script_dir + "/lib/RNA/ACGU.rna.pdb", check_dna=True)
    pairs = [[0, 7], [1, 6], [2, 5], [3, 4]]
    forms = helix_form(query_atom23_pos[:, [0, 5, 12], :], pairs, script_dir)
    print(forms)

    query_atom23_pos, _, _, _, _, _ = read_atom23_na(script_dir + "/lib/DNA/ACGT.dna.pdb", check_dna=True)
    pairs = [[0, 7], [1, 6], [2, 5], [3, 4]]
    forms = helix_form(query_atom23_pos[:, [0, 5, 12], :], pairs, script_dir)
    print(forms)


