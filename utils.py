import os
import sys
import mrcfile
import warnings
import numpy as np
from math import floor, ceil
try:
    from interp3d import interp3d as f90interp3d
    can_import_f90interp3d = True
except Exception as e:
    can_import_f90interp3d = False
    raise e

'''parse_map, pad_map, split_map_into_overlapped_chunks, get_map_from_overlapped_chunks, write_map'''
def split_map_into_overlapped_chunks(map, box_size, stride, dtype=np.float32, padding=0.0):
    assert stride < box_size
    map_shape = np.shape(map)
    padded_map = np.full((map_shape[0] + 2 * box_size, map_shape[1] + 2 * box_size, map_shape[2] + 2 * box_size), padding, dtype=dtype)
    padded_map[box_size : box_size + map_shape[0], box_size : box_size + map_shape[1], box_size : box_size + map_shape[2]] = map
    chunk_list = list()
    start_point = box_size - stride
    cur_x, cur_y, cur_z = start_point, start_point, start_point
    while (cur_z + stride < map_shape[2] + box_size):
        next_chunk = padded_map[cur_x:cur_x + box_size, cur_y:cur_y + box_size, cur_z:cur_z + box_size]
        cur_x += stride
        if (cur_x + stride >= map_shape[0] + box_size):
            cur_y += stride
            cur_x = start_point # Reset
            if (cur_y + stride  >= map_shape[1] + box_size):
                cur_z += stride
                cur_y = start_point # Reset
                cur_x = start_point # Reset
        chunk_list.append(next_chunk)
    n_chunks = len(chunk_list)
    ncx, ncy, ncz = [ceil(map_shape[i] / stride) for i in range(3)]
    assert(n_chunks == ncx * ncy * ncz)
    chunks = np.asarray(chunk_list, dtype=dtype)
    return chunks, ncx, ncy, ncz

def pad_map(map, box_size, core_size, dtype=np.float32, padding=0):
    map_shape = np.shape(map)
    ncx, ncy, ncz = [ceil(map_shape[i] / core_size) for i in range(3)]

    padded_map = np.full((ncx * core_size, ncy * core_size, ncz * core_size), padding, dtype=dtype)
    padded_map[:map_shape[0], :map_shape[1], :map_shape[2]] = map
    
    return padded_map
    
def get_map_from_overlapped_chunks(chunks, ncx, ncy, ncz, n_classes, box_size, stride, nxyz, dtype=np.float32):
    map = np.zeros((n_classes, \
                    (ncx - 1) * stride + box_size, \
                    (ncy - 1) * stride + box_size, \
                    (ncz - 1) * stride + box_size), dtype=dtype)
    denominator = np.zeros((n_classes, \
                            (ncx - 1) * stride + box_size, \
                            (ncy - 1) * stride + box_size, \
                            (ncz - 1) * stride + box_size), dtype=dtype) # should clip to 1
    i = 0
    for z_steps in range(ncz):
        for y_steps in range(ncy):
            for x_steps in range(ncx):
                map[:, x_steps * stride : x_steps * stride + box_size,
                       y_steps * stride : y_steps * stride + box_size,
                       z_steps * stride : z_steps * stride + box_size] += chunks[i]
                denominator[:, x_steps * stride : x_steps * stride + box_size,
                               y_steps * stride : y_steps * stride + box_size,
                               z_steps * stride : z_steps * stride + box_size] += 1

                i += 1
    return (map / denominator.clip(min=1))[:, stride : nxyz[2] + stride, stride : nxyz[1] + stride, stride : nxyz[0] + stride]


def parse_map(map_file, ignorestart, apix=None, origin_shift=None, interp='f90'):
    assert interp in ['py', 'f90', 'fortran']
    print("# Intend to use interp package -> '{}'".format(interp))
    if interp in ['f90', 'fortran']:
        interp3d = f90interp3d
    else:
        from pyinterp3d import pyinterp3d
        interp3d = pyinterp3d

    ''' parse mrc '''
    mrc = mrcfile.open(map_file, mode='r')

    map = np.asfarray(mrc.data.copy(), dtype=np.float32)
    voxel_size = np.asarray([mrc.voxel_size.x, mrc.voxel_size.y, mrc.voxel_size.z], dtype=np.float32)
    ncrsstart = np.asarray([mrc.header.nxstart, mrc.header.nystart, mrc.header.nzstart], dtype=np.float32)
    origin = np.asarray([mrc.header.origin.x, mrc.header.origin.y, mrc.header.origin.z], dtype=np.float32)
    ncrs = (mrc.header.nx, mrc.header.ny, mrc.header.nz)
    angle = np.asarray([mrc.header.cellb.alpha, mrc.header.cellb.beta, mrc.header.cellb.gamma], dtype=np.float32)

    ''' check orthogonal '''
    try:
        assert(angle[0] == angle[1] == angle[2] == 90.0)
    except AssertionError:
        print("# Input grid is not orthogonal. EXIT.")
        mrc.close()
        sys.exit()

    ''' reorder axes '''
    mapcrs = np.subtract([mrc.header.mapc, mrc.header.mapr, mrc.header.maps], 1)
    sort = np.asarray([0, 1, 2], dtype=np.int32)
    for i in range(3):
        sort[mapcrs[i]] = i
    nxyzstart = np.asarray([ncrsstart[i] for i in sort])
    nxyz = np.asarray([ncrs[i] for i in sort])

    map = np.transpose(map, axes=2-sort[::-1])
    mrc.close()

    ''' shift origin according to n*start '''
    if not ignorestart:
        origin += np.multiply(nxyzstart, voxel_size)

    ''' shift by decimal '''
    if origin_shift is not None: 
        origin_shift = origin_shift - origin + np.floor(origin)

    ''' interpolate grid interval '''
    if apix is not None:
        try:
            assert(voxel_size[0] == voxel_size[1] == voxel_size[2] == apix and origin_shift is None)
        except AssertionError:
            interp3d.del_mapout()
            target_voxel_size = np.asarray([apix, apix, apix], dtype=np.float32)
            print("# Rescale voxel size from {} to {}, shift origin by {}".format(voxel_size, target_voxel_size, origin_shift))
            if origin_shift is not None:
                interp3d.cubic(map, voxel_size[2], voxel_size[1], voxel_size[0], apix, origin_shift[2], origin_shift[1], origin_shift[0], nxyz[2], nxyz[1], nxyz[0])
                origin += origin_shift
            else:
                interp3d.cubic(map, voxel_size[2], voxel_size[1], voxel_size[0], apix, 0.0, 0.0, 0.0, nxyz[2], nxyz[1], nxyz[0])
                
            map = interp3d.mapout
            nxyz = np.asarray([interp3d.pextx, interp3d.pexty, interp3d.pextz], dtype=np.int32)
            voxel_size = target_voxel_size

    assert(np.all(nxyz == np.asarray([map.shape[2], map.shape[1], map.shape[0]], dtype=np.int32)))

    return map, origin, nxyz, voxel_size

def write_map(file_name, map, voxel_size, origin=(0.0, 0.0, 0.0), nxyzstart=(0, 0, 0)):
    mrc = mrcfile.new(file_name, overwrite=True)
    mrc.set_data(map)
    (mrc.header.nxstart, mrc.header.nystart, mrc.header.nzstart) = nxyzstart
    (mrc.header.origin.x, mrc.header.origin.y, mrc.header.origin.z) = origin
    mrc.voxel_size = [voxel_size[i] for i in range(3)]

    mrc.close()

if __name__ == "__main__":
    import sys
    data0, origin0, _, vsize0 = parse_map(sys.argv[1], False, 1.0, interp='f90')
    write_map("f90.mrc", data0, vsize0, origin0)
    data1, origin1, _, vsize1 = parse_map(sys.argv[1], False, 1.0, interp='py')
    write_map("py.mrc",  data1, vsize1, origin1)
    deviation = np.abs(data0-data1).mean()
    print("# Overall floating point deviation = {}".format(deviation))

