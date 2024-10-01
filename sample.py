import os
import sys
import math
import time
import argparse
import tempfile
import numpy as np
from utils import parse_map
from scipy.spatial import KDTree
from tsptrace import trace
from resample import upsample_3d_points, smooth
from pdbio import read_pdb, flatten
# Set random seed
np.random.seed(0)

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=str)
    parser.add_argument("--c4", type=str)
    parser.add_argument("--mc", type=str)

    parser.add_argument("--pmap", type=str)
    parser.add_argument("--c4map", type=str)

    parser.add_argument("--lkhdir", type=str)
    parser.add_argument("--outdir", type=str)

    parser.add_argument("--length", "-l", type=int)

    args = parser.parse_args()
    return args

def write_points_as_pdb(f, coords, aatypes=None, bfactors=None, ter=False, end=False, model=1, chain='A', atomtype=' CA ', mode='w'):
    if len(coords) == 0:
        return

    if bfactors is None or len(coords) != len(bfactors):
        bfactors = [0.0 for _ in range(len(coords))]
    
    if aatypes is None or len(coords) != len(aatypes):
        aatypes = ['U' for _ in range(len(coords))]

    original_stdout = sys.stdout

    if model is not None:
        f.write(f"MODEL{model:>9d}\n")

    for i, coord in enumerate(coords):
        f.write("ATOM  {:5d} {:4s}   {:1s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n".format(i+1, atomtype, aatypes[i], chain, i+1, coord[0], coord[1], coord[2], 1.0, bfactors[i]))
        if i < len(coords) and ter:
            fwrite("TER\n")
    if not ter:
        f.write("TER\n")

    if end:
        f.write("END\n")


def safe_norm(v, eps=1e-3):
    v += 1e-3
    norm = max(eps, np.linalg.norm(v))
    return v / norm

def get_density(grid, coord):
    x, y, z = [int(x) for x in coord]
    n0, n1, n2 = grid.shape
    if x >= n2 or x < 0 or \
       y >= n1 or y < 0 or \
       z >= n0 or z < 0:
        return 0.
    else:
        return grid[z, y, x]

def distance(x, y):
    return np.sqrt(np.sum(np.power(np.subtract(x, y), 2)))

def filter_with_distance(chain, density, r0=3.):
    while True:
        tree = KDTree(chain)
        keep = [False for _ in range(len(chain))]
        for i in range(len(chain)):
            inds = tree.query_ball_point(chain[i], r0)
            dens = density[inds]

            ii = np.argmax(dens)
            ii = inds[ii]
            keep[ii] = True
        if np.all(keep):
            break

        chain = chain[keep]
        density = density[keep]

    return chain

def is_near_end(index, size, l=4):
    return index <= l or abs(size - index) <= l

def get_subarrays(arr):
    subarrays = []
    subarray = [arr[0]]

    for i in range(1, len(arr)):
        if arr[i] == arr[i - 1] + 1:
            subarray.append(arr[i])
        else:
            subarrays.append(subarray)
            subarray = [arr[i]]
    
    subarrays.append(subarray)
    return subarrays

def filter_paths(paths):
    used_paths = set()
    single_paths = []
    for i, p in enumerate(paths):
        # We consider forward/backward as same path
        p_forward = tuple(p)
        p_backward = tuple(p[::-1])
        if (p_forward in used_paths) or (p_backward in used_paths):
            continue
        used_paths.add(p_forward)
        used_paths.add(p_backward)
        single_paths.append(p)
    paths = single_paths
    return paths


def split_chain_by_breaks(chain, breaks):
    ends = list(breaks)
    ends.append(len(chain)-1)
    ends.sort()
    inds = [i for i in range(len(chain))]
    start = 0
    segs = []
    for end in ends:
        segs.append(inds[start:end+1])
        start = end+1
    if len(segs) == 0:
        segs = [inds]
    chains = []
    for seg in segs:
        chains.append(chain[seg])
    return chains


def equal_sample_and_infill(chain):
    last = 0
    dc4 = 5.
    d0 = 0.
    d1 = 5.
    pruned_chain = [last]
    while last + 1 < len(chain):
        next = last + 1
        while next < len(chain):
            d = distance(chain[last], chain[next])
            print(d)
            if (d0 < d < d1):
                next += 1
                continue
            else:
                break
        next0 = next - 1
        if next0 == last:
            # We should add points
            d = distance(last, last + 1)
            n = int(d / dc4)
            pruned_chain.append(next0 + 1)
            last = next0 + 1
            for k in range(n):
                pass
        else:
            # We should use next0 as last's next node
            print(f"Add {next0} to {last}")
            pruned_chain.append(next0)
            last = next0
    return pruned_chain

