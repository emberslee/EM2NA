import os
import numpy as np
from numba import float32, int32, float64
from numba import jit

# Max thread
os.environ['NUMBA_NUM_THREADS'] = '4'

# For pyinterp3d.cubic
@jit(nopython=True)
def get_w(x):
    a = -0.5
    w = np.zeros(4)
    intx = np.floor(x)
    d1 = 1.0 + (x - intx)
    d2 = d1 - 1.0
    d3 = 1.0 - d2
    d4 = d3 + 1.0
    
    w[0] = a * np.abs(d1**3) - 5 * a * d1**2 + 8 * a * np.abs(d1) - 4 * a
    w[1] = (a + 2) * np.abs(d2**3) - (a + 3) * d2**2 + 1
    w[2] = (a + 2) * np.abs(d3**3) - (a + 3) * d3**2 + 1
    w[3] = a * np.abs(d4**3) - 5 * a * d4**2 + 8 * a * np.abs(d4) - 4 * a
    return w

class pyinterp3d:
    pextx = None
    pexty = None
    pextz = None
    mapout = None

    @staticmethod
    def del_mapout():
        pextx = None
        pexty = None
        pextz = None
        mapout = None

    @staticmethod
    @jit(nopython=True, parallel=True)
    def _linear(mapin, zpix, ypix, xpix, apix, shiftz, shifty, shiftx, nz, ny, nx):
        pextx0 = int(np.floor(xpix * (nx - 1) / apix)) + 1
        pexty0 = int(np.floor(ypix * (ny - 1) / apix)) + 1
        pextz0 = int(np.floor(zpix * (nz - 1) / apix)) + 1
        mapout0 = np.zeros((pextz0, pexty0, pextx0), dtype=mapin.dtype)
        
        for indz in range(pextz0):
            for indy in range(pexty0):
                for indx in range(pextx0):
                    xpos = indx * apix + shiftx
                    ypos = indy * apix + shifty
                    zpos = indz * apix + shiftz
                    gx = xpos / xpix + 1
                    gy = ypos / ypix + 1
                    gz = zpos / zpix + 1
                    x0, y0, z0 = int(np.floor(gx)) - 1, int(np.floor(gy)) - 1, int(np.floor(gz)) - 1
                    x1, y1, z1 = x0 + 1, y0 + 1, z0 + 1
                    if x0 >= 0 and x1 < nx and y0 >= 0 and y1 < ny and z0 >= 0 and z1 < nz:
                        a, b, c = gx - (x0 + 1), gy - (y0 + 1), gz - (z0 + 1)
                        mapout0[indz, indy, indx] = (
                            a * b * c * mapin[z1, y1, x1] +
                            (1-a) * b * c * mapin[z1, y1, x0] +
                            a * (1-b) * c * mapin[z1, y0, x1] +
                            a * b * (1-c) * mapin[z0, y1, x1] +
                            a * (1-b) * (1-c) * mapin[z0, y0, x1] +
                            (1-a) * b * (1-c) * mapin[z0, y1, x0] +
                            (1-a) * (1-b) * c * mapin[z1, y0, x0] +
                            (1-a) * (1-b) * (1-c) * mapin[z0, y0, x0]
                        )

        return mapout0, pextx0, pexty0, pextz0

    @staticmethod
    @jit(nopython=True, parallel=True)
    def _cubic(mapin, zpix, ypix, xpix, apix, shiftz, shifty, shiftx, nz, ny, nx):
        pextx0 = int(np.floor(xpix * (nx - 1) / apix)) + 1
        pexty0 = int(np.floor(ypix * (ny - 1) / apix)) + 1
        pextz0 = int(np.floor(zpix * (nz - 1) / apix)) + 1
        
        mapout0 = np.zeros((pextz0, pexty0, pextx0), dtype=mapin.dtype)
        
        for indz in range(pextz0):
            for indy in range(pexty0):
                for indx in range(pextx0):
                    gz = (indz * apix + shiftz) / zpix
                    gy = (indy * apix + shifty) / ypix
                    gx = (indx * apix + shiftx) / xpix
                    
                    intz = int(np.floor(gz))
                    inty = int(np.floor(gy))
                    intx = int(np.floor(gx))
                    
                    wz = get_w(gz - intz)
                    wy = get_w(gy - inty)
                    wx = get_w(gx - intx)
                    
                    for i in range(4):
                        for j in range(4):
                            for k in range(4):
                                zi = intz + i - 1
                                yi = inty + j - 1
                                xi = intx + k - 1
                                if 0 <= zi < nz and 0 <= yi < ny and 0 <= xi < nx:
                                    mapout0[indz, indy, indx] += mapin[zi, yi, xi] * wz[i] * wy[j] * wx[k]

        return mapout0, pextx0, pexty0, pextz0

    # Interfaces
    @staticmethod
    def linear(mapin, zpix, ypix, xpix, apix, shiftz, shifty, shiftx, nz, ny, nx):
        mapout0, pextx0, pexty0, pextz0 = pyinterp3d._linear(mapin, zpix, ypix, xpix, apix, shiftz, shifty, shiftx, nz, ny, nx)
        pyinterp3d.mapout = mapout0
        pyinterp3d.pextx = pextx0
        pyinterp3d.pexty = pexty0
        pyinterp3d.pextz = pextz0
        return mapout0

    @staticmethod
    def cubic(mapin, zpix, ypix, xpix, apix, shiftz, shifty, shiftx, nz, ny, nx):
        mapout0, pextx0, pexty0, pextz0 = pyinterp3d._cubic(mapin, zpix, ypix, xpix, apix, shiftz, shifty, shiftx, nz, ny, nx)
        pyinterp3d.mapout = mapout0
        pyinterp3d.pextx = pextx0
        pyinterp3d.pexty = pexty0
        pyinterp3d.pextz = pextz0
        return mapout0


