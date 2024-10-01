import numpy as np

def get_aa_type(grid : np.ndarray, coord : np.ndarray, origin=None):
    # AGCU/T
    if origin is not None:
        coord = np.subtract(coord, origin)
    x, y, z = [int(x) for x in coord]
    n0, n1, n2 = grid.shape
    if x >= n2 or x < 0 or \
       y >= n1 or y < 0 or \
       z >= n0 or z < 0:
        return 3
    else:
        # Vote for 3**3 voxels
        votes = [0] * 4
        for xx in range(x-1, x+2):
            for yy in range(y-1, y+2):
                for zz in range(z-1, z+2):
                    if 0 <= xx < n2 and 0 <= yy < n1 and 0 <= zz < n0:
                        v = int(grid[zz, yy, xx])
                        if v < 4:
                            votes[v] += 1
        return np.argmax(votes)


def get_grid_value(grid : np.ndarray, coord : np.ndarray, origin=None):
    if origin is not None:
        coord = np.subtract(coord, origin)
    x, y, z = [int(x) for x in coord]
    n0, n1, n2 = grid.shape
    if x >= n2 or x < 0 or \
       y >= n1 or y < 0 or \
       z >= n0 or z < 0:
        return 0.0
    else:
        return grid[z, y, x]


def get_grid_value_interp(grid : np.ndarray, coord : np.ndarray, origin=None):
    n0, n1, n2 = grid.shape

    if origin is not None:
        coord = np.subtract(coord, origin)

    z = 0.0

    kx0 = int(coord[0])
    kx1 = kx0 + 1

    ky0 = int(coord[1])
    ky1 = ky0 + 1

    kz0 = int(coord[2])
    kz1 = kz0 + 1

    if  kx0 < 0 or kx0 >= n2 or \
        kx1 < 0 or kx1 >= n2 or \
        ky0 < 0 or ky0 >= n1 or \
        ky1 < 0 or ky1 >= n1 or \
        kz0 < 0 or kz0 >= n0 or \
        kz1 < 0 or kz1 >= n0:
        return 0.0

    x0 = coord[0] - kx0
    y0 = coord[1] - ky0
    z0 = coord[2] - kz0

    txy = x0*y0
    tyz = y0*z0
    txz = x0*z0
    txyz = x0*y0*z0

    # z, y, x
    v000 = grid[kz0, ky0, kx0]
    v100 = grid[kz1, ky0, kx0]
    v010 = grid[kz0, ky1, kx0]
    v001 = grid[kz0, ky0, kx1]
    v101 = grid[kz1, ky0, kx1]
    v011 = grid[kz0, ky1, kx1]
    v110 = grid[kz1, ky1, kx0]
    v111 = grid[kz1, ky1, kx1]

    temp1 = v000*(1.0-x0-y0-z0+txy+tyz+txz-txyz);
    temp2 = v100*(x0-txy-txz+txyz);
    temp3 = v010*(y0-txy-tyz+txyz);
    temp4 = v001*(z0-txz-tyz+txyz);
    temp5 = v101*(txz-txyz);
    temp6 = v011*(tyz-txyz);
    temp7 = v110*(txy-txyz);
    temp8 = v111*txyz;

    z = temp1 + temp2 + temp3 + temp4 + temp5 + temp6 + temp7 + temp8

    return z
