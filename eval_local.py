import sys
import numpy as np
import argparse

from grid import (
    get_grid_value,
    get_grid_value_interp,
)

from pdbio import (
    nuc_type_index_to_1,
    nuc_type_1_to_index,
    read_atom23_na,
    res_chains_atom23_to_pdb,
)

from utils import (
    parse_map,
)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", "-p", type=str, required=True, help="Input structure")
    parser.add_argument("--out", "-o", type=str, default="./score.pdb", help="Output structure")
    parser.add_argument("--pmap", "-pmap", type=str, required=True, help="Predicted P map")
    parser.add_argument("--cmap", "-cmap", type=str, required=True, help="Predicted C map")
    parser.add_argument("--nmap", "-nmap", type=str, required=True, help="Predicted N map")
    parser.add_argument("--atom", "-a", type=str, default="C4", help="Atom selection, P or C4' or N")
    parser.add_argument("--percentile", "-percentile", type=float, default=99.9, help="Map percentile as maximum to normalize score")
    parser.add_argument("--nointerp", action='store_true', help="Not use interpolated density")
    # Optional interp
    parser.add_argument("--pyinterp", action='store_true', help="Interpolation using pyinterp3d")
    args = parser.parse_args()
    return args


def smooth_window(data, kernel=[1.0, 4.0, 9.0, 4.0, 1.0]):
    # Kernel size must be odd
    assert len(kernel) % 2 == 1, "# Error Kernel size must be odd number"

    kernel = np.asarray(kernel, dtype=np.float32)
    window = len(kernel) // 2

    # . . . . . . . c . . . . . . .
    # . . . . . w w w w w . . . . .
    ret = np.zeros(len(data), dtype=np.float32)
    L = len(data)
    for center in range(L):
        l = center - window
        r = center + window + 1
        s = 0
        w = 1e-3

        for k in range(l, r):
            if 0 <= k < len(data):
                s += data[k] * kernel[k - l]
                w += kernel[k - l]
        ret[center] = s / w
    return ret

if __name__ == '__main__':
    args = get_args()
    fout = args.out
    fpdb = args.pdb
    fpmap = args.pmap
    fcmap = args.cmap
    fnmap = args.nmap
    dpercentile = args.percentile
    if not (99.0 <= dpercentile <= 100.0):
        print("# WARN Invalid percentile must be in [99.0, 100.0]")
        dpercentile = 99.9
        print("# WARN Enforce percentile to be {:.4f}".format(dpercentile))
    else:
        print("# Use a percentile of {:.4f} to normalize the density".format(dpercentile))


    # Read map
    if args.pyinterp:
        interp = "py"
    else:
        interp = "f90"
    pmap, porigin, pnxyz, pvsize = parse_map(fpmap, False, None, interp=interp)
    cmap, corigin, cnxyz, cvsize = parse_map(fcmap, False, None, interp=interp)
    nmap, norigin, nnxyz, nvsize = parse_map(fnmap, False, None, interp=interp)

    atom_select = args.atom
    atom_index = 5
    sgrid = cmap
    origin = corigin
    vsize = cvsize
    if atom_select in ["C4", "C4'", "C", "c4", "c4'", "c"]:
        atom_index = 5
        sgrid = cmap
        origin = corigin
        vsize = cvsize
    elif atom_select in ["P", "p"]:
        atom_index = 0
        sgrid = pmap
        origin = porigin
        vsize = pvsize
    elif atom_select in ["N", "N1", "N9", "n", "n1", "n9"]:
        atom_index = 12
        sgrid = nmap
        origin = norigin
        vsize = nvsize
    else:
        print("# WARN Wrong atom selection use C4' instead")
        atom_index = 5
        sgrid = cmap
        origin = corigin
        vsize = nvsize

    # Read PDB
    atom23_pos, atom23_mask, res_types, res_indices, chain_indices, res_is_dna = read_atom23_na(fpdb, check_dna=True)


    # Normalize sgrid
    smin = 0.0
    smax = np.percentile(sgrid[sgrid > 0.0], q=dpercentile)
    sgrid = np.clip(sgrid, a_min=smin, a_max=smax)
    sgrid = sgrid / smax

    # Get score
    atom_pos = atom23_pos[:, atom_index, :]
    get_grid_value_func = get_grid_value_interp
    if args.nointerp:
        get_grid_value_func = get_grid_value

    dens = [get_grid_value_func(sgrid, (atom_pos[i]-origin)/vsize) for i in range(len(atom_pos))]
    #print(dens)

    # Use dens as bfactors
    bfactors = np.asarray(dens, dtype=np.float32)
    bfactors = np.sqrt(bfactors) * (bfactors >= 0.5) + bfactors * (bfactors < 0.5)


    #bfactors = np.clip(bfactors, a_min=0.00, a_max=0.99)
    assert len(bfactors) == len(atom23_pos)


    chains_atom23_pos = []
    chains_atom23_mask = []
    chains_res_types = []
    chains_bfactors = []
    chains_is_dna = []

    n_chains = np.max(chain_indices) + 1
    all_smoothed_bfactors = []

    for i in range(0, n_chains):
        chains_atom23_pos.append(atom23_pos[chain_indices == i])
        chains_atom23_mask.append(atom23_mask[chain_indices == i])
        chains_res_types.append([nuc_type_index_to_1[x] for x in res_types[chain_indices == i]])
        raw_bfactors = bfactors[chain_indices == i]
        # Smooth bfactor using shifted window
        smoothed_bfactors = smooth_window(raw_bfactors)
        smoothed_bfactors = np.clip(smoothed_bfactors, a_min=0.00, a_max=0.99)

        # 0702
        smoothed_bfactors = (smoothed_bfactors * 10000).astype(np.int32) / 10000

        all_smoothed_bfactors.append(smoothed_bfactors)

        chains_bfactors.append(smoothed_bfactors)
        chains_is_dna.append(res_is_dna[chain_indices == i])


    all_smoothed_bfactors = np.concatenate(all_smoothed_bfactors, axis=0)

    avg_bfactor = np.mean(all_smoothed_bfactors)
    med_bfactor = np.median(all_smoothed_bfactors)
    max_bfactor = np.max(all_smoothed_bfactors)
    min_bfactor = np.min(all_smoothed_bfactors)
    print("# Report local quality estimation")
    print("# Median  bfactor/plddt is {:.4f}".format(med_bfactor))
    print("# Average bfactor/plddt is {:.4f}".format(avg_bfactor))
    print("# Maximum bfactor/plddt is {:.4f}".format(max_bfactor))
    print("# Minimum bfactor/plddt is {:.4f}".format(min_bfactor))

    res_chains_atom23_to_pdb(
        filename=fout,
        chains_atom23_pos=chains_atom23_pos,
        chains_atom23_mask=chains_atom23_mask,
        chains_res_types=chains_res_types,
        chains_occupancy=None,
        chains_bfactors=chains_bfactors,
        chains_is_dna=chains_is_dna,
    )
    print("# Write scored PDB to {}".format(fout))


