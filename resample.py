import math
import argparse
import numpy as np
from scipy.spatial import KDTree
from scipy.interpolate import CubicSpline
from pdbio import flatten, read_pdb

def upsample_3d_points(points, upsample_factor):
    num_points = len(points)
    t = np.linspace(0, 1, num_points)
    t_upsampled = np.linspace(0, 1, int(num_points * upsample_factor))

    # Separate x, y, z coordinates of the points
    x = [p[0] for p in points]
    y = [p[1] for p in points]
    z = [p[2] for p in points]

    # Create cubic spline interpolators for x, y, and z coordinates
    spline_x = CubicSpline(t, x)
    spline_y = CubicSpline(t, y)
    spline_z = CubicSpline(t, z)

    # Interpolate the upsampled points
    x_upsampled = spline_x(t_upsampled)
    y_upsampled = spline_y(t_upsampled)
    z_upsampled = spline_z(t_upsampled)

    # Combine the upsampled coordinates with the original coordinates to form the new point cloud
    combined_points = []
    for i in range(num_points - 1):
        combined_points.append((x[i], y[i], z[i]))
        j = 1
        while j < upsample_factor:
            t_sample = t[i] + j * (t[i + 1] - t[i]) / upsample_factor
            combined_points.append((spline_x(t_sample), spline_y(t_sample), spline_z(t_sample)))
            j += 1

    return combined_points


def smooth(chain):
    chain_copy = [p for p in chain]
    n_round = 0
    has_sharp_frag = True
    while True:
       if not has_sharp_frag:
           break
       has_sharp_frag = False

       chain0 = [chain[0]]
       i = 1
       while i + 1 < len(chain):
           # Consider angles between i - 1, i, i + 1
           k = i
           while k + 1 < len(chain):
               v1 = chain[i - 1] - chain[k]
               v2 = chain[k + 1] - chain[k]
               sign = np.inner(v1, v2)
               if sign > 0:
                   has_sharp_frag = True
                   k += 1
               else:
                   break
           # Skip from [i, k]
           chain0.append(chain[k])
           i = k + 1

       chain0.append(chain[-1])
       # Update
       chain = np.asarray(chain0)
       n_round += 1
       print(f"Smooth round {n_round}")

   # Split chain based on probability
    chain = chain.tolist()
    return chain

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mc", "-mc", type=str)
    parser.add_argument("--p1", "-p1", type=str)
    parser.add_argument("--p2", "-p2", type=str)
    args = parser.parse_args()
    return args

def print_points_as_pdb(coords, bfactors=None, ter=False, end=False, model=1):
    if len(coords) == 0:
        return
    if bfactors is None or len(coords) != len(bfactors):
        bfactors = [0.0 for _ in range(len(coords))]
    
    if model is not None:
        print(f"MODEL{model:>9d}")

    for i, coord in enumerate(coords):
        print("ATOM  {:5d}  CA    U A{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}".format(i+1, i+1, coord[0], coord[1], coord[2], 1.0, bfactors[i]))
        if i < len(coords) and ter:
            print("TER")
    if not ter:
        print("TER")    

    if end:
        print("END")


def safe_norm(v, eps=1e-3):
    v += 1e-3
    norm = max(eps, np.linalg.norm(v))
    return v / norm


