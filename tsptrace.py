import os
import tempfile
import numpy as np
#from sklearn.metrics import pairwise_distances
from geo import pairwise_distances

def set_tsp_file(mat, path='./', tsptype='TSP'):
    # Write to file
    filename = path + '/exp.tsp'
    f = open(filename, 'w')
    f.write('NAME : thread\nCOMMENT : thread\nTYPE : {}\nDIMENSION : {}\nEDGE_WEIGHT_TYPE : EXPLICIT\nEDGE_WEIGHT_FORMAT : FULL_MATRIX\nEDGE_WEIGHT_SECTION\n'.format(tsptype, len(mat)))

    for i in range(len(mat) - 1):
        for j in range(len(mat) - 1):
            f.write('{} '.format(int(mat[i][j])))
        f.write('0\n')
    for i in range(len(mat)):
        f.write('0 ')
    f.write('\nEOF\n')
    #print('Write to tsp file to {}'.format(path))


def set_tsp_params(depot=1, vehicles=1, path='./', time_limit=None):
    # Write params to file
    filename = path + '/exp.par'
    with open(filename, 'w') as f:
        if isinstance(time_limit, int) and time_limit > 0:
            f.write('PROBLEM_FILE = {}\nOPTIMUM = 99999999\nMOVE_TYPE = 5\nPATCHING_C = 3\nPATCHING_A = 2\nRUNS = 5\nVEHICLES = {}\nDEPOT = {}\nTOUR_FILE = {}\nTIME_LIMIT = {}\n'.format( os.path.join(path, 'exp.tsp'), vehicles, depot, os.path.join(path, 'exp.out'), time_limit))
        else:
            f.write('PROBLEM_FILE = {}\nOPTIMUM = 99999999\nMOVE_TYPE = 5\nPATCHING_C = 3\nPATCHING_A = 2\nRUNS = 5\nVEHICLES = {}\nDEPOT = {}\nTOUR_FILE = {}\n'.format( os.path.join(path, 'exp.tsp'), vehicles, depot, os.path.join(path, 'exp.out')))
    #print('Write tsp params to {}'.format(path))


def run_tsp(vehicles=1, path='./', lkh_dir='./'):
    filename = path + '/exp.par'
    cmd = lkh_dir + '/LKH ' + filename + ' > /dev/null'
    os.system(cmd)
    #print('Run LKH for num path = {}'.format(vehicles))


def distance(x, y):
    return np.sqrt(np.sum(np.power(np.subtract(x, y), 2)))

def trace(points, vehicles=1, lkh='./', time_limit=None):
    # Get abs dir for lkh
    lkh_dir = os.path.abspath(  os.path.expanduser(lkh)  )
    #print("Set lkh path to : ", lkh_dir)
    scale_factor = 10.
    # Make a temperary dir
    with tempfile.TemporaryDirectory() as temp_dir:
        #print("Temporary directory created : ", temp_dir)

        # Set lkh parameters
        L = len(points) + 1
        dmat = np.zeros((L, L), dtype=np.float32)
        dmat[:L-1, :L-1] = pairwise_distances(points, points) * scale_factor

        set_tsp_file(dmat, path=temp_dir)
        set_tsp_params(depot=len(points) + 1, vehicles=vehicles, path=temp_dir, time_limit=time_limit)
        #print("Set up tsp files")

        # Trace on lkh
        run_tsp(vehicles=vehicles, path=temp_dir, lkh_dir=lkh_dir)
        print("# Running lkh")

        # Get result
        with open(os.path.join(temp_dir, 'exp.out'), 'r') as f:
            lines = f.readlines()
        sol = list()
        flag = False
        for line in lines:
            if line.startswith('TOUR_SECTION'):
                flag = True
                continue
            if line.startswith('-1'):
                break
            if flag:
                x = int(line.strip().split()[0])
                sol.append(x)
        out_inds = list()
        all_out_inds = list()
        for i, ind in enumerate(sol):
            if 1 <= ind <= len(points):
                out_inds.append(ind - 1)
            else:
                all_out_inds.append(out_inds)
                out_inds = []
        all_out_inds.append(out_inds)
        all_out_inds = all_out_inds[1:]

    # Return traced points (indices)
    return all_out_inds


def residue_trace(residues, vehicles=1, lkh='./'):
    # Get abs dir for lkh
    lkh_dir = os.path.abspath(  os.path.expanduser(lkh)  )
    scale_factor = 10.
    L = len(residues) + 1
    n_positions = residues[:, 0, :]
    c_positions = residues[:, 2, :]
    dmat = np.zeros((L, L), dtype=np.float32)
    dmat[:L-1, :L-1] = pairwise_distances(n_positions, c_positions) * scale_factor
    # Make a temperary dir
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = './'
        #print("Temporary directory created : ", temp_dir)

        # Set lkh parameters
        set_tsp_file(dmat, path=temp_dir, tsptype='ATSP')
        set_tsp_params(depot=L, vehicles=vehicles, path=temp_dir)
        #print("Set up tsp files")

        # Trace on lkh
        run_tsp(vehicles=vehicles, path=temp_dir, lkh_dir=lkh_dir)
        print("Running fragment lkh")

        # Get result
        with open(os.path.join(temp_dir, 'exp.out'), 'r') as f:
            lines = f.readlines()

        sol = list()
        flag = False
        for line in lines:
            if line.startswith('TOUR_SECTION'):
                flag = True
                continue
            if line.startswith('-1'):
                break
            if flag:
                x = int(line.strip().split()[0])
                sol.append(x)
        out_inds = list()
        all_out_inds = list()
        for i, ind in enumerate(sol):
            if 1 <= ind <= L - 1:
                out_inds.append(ind - 1)
            else:
                all_out_inds.append(out_inds)
                out_inds = []
        all_out_inds.append(out_inds)
        all_out_inds = all_out_inds[1:]

    # Return traced
    return all_out_inds
    


# Trace fragment
def fragment_trace(
        frags, 
        vehicles=1,
        lkh='./',
        ordered=True
    ):
    def distance(x, y):
        return np.sqrt(np.sum(np.power(np.subtract(x, y), 2)))

    # Special cases
    if len(frags) == 1:
        return [[0]]
    if len(frags) == 2:
        d0 = distance(frags[0][-1][-1], frags[-1][0][0])
        d1 = distance(frags[-1][-1][-1], frags[0][0][0])
        if d0 < d1:
            return [[0, 1]]
        else:
            return [[1, 0]]

    # Get abs dir for lkh
    lkh_dir = os.path.abspath(  os.path.expanduser(lkh)  )
    #print("Set lkh path to : ", lkh_dir)

    # True distance between fragment terminals
    # Pseudo distance between 2 terminal in a fragment
    npoint = len(frags)
    dmat = np.zeros((npoint + 1, npoint + 1), dtype=np.float32)
    scale_factor = 10.

    # Only start from head to foot
    for i in range(len(frags)):
        fi = frags[i]
        for j in range(len(frags)):
            fj = frags[j]
            # i to j
            dmat[i, j] = distance(fi[-1][-1], fj[0][0]) * scale_factor

    # Make a temperary dir
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = './'
        #print("Temporary directory created : ", temp_dir)

        # Set lkh parameters
        set_tsp_file(dmat, path=temp_dir, tsptype='ATSP')
        set_tsp_params(depot=npoint + 1, vehicles=vehicles, path=temp_dir)
        #print("Set up tsp files")

        # Trace on lkh
        run_tsp(vehicles=vehicles, path=temp_dir, lkh_dir=lkh_dir)
        print("Running fragment lkh")

        # Get result
        with open(os.path.join(temp_dir, 'exp.out'), 'r') as f:
            lines = f.readlines()

        sol = list()
        flag = False
        for line in lines:
            if line.startswith('TOUR_SECTION'):
                flag = True
                continue
            if line.startswith('-1'):
                break
            if flag:
                x = int(line.strip().split()[0])
                sol.append(x)
        out_inds = list()
        all_out_inds = list()
        for i, ind in enumerate(sol):
            if 1 <= ind <= npoint:
                out_inds.append(ind - 1)
            else:
                all_out_inds.append(out_inds)
                out_inds = []
        all_out_inds.append(out_inds)
        all_out_inds = all_out_inds[1:]

    # Return traced
    return all_out_inds


if __name__ == '__main__':
    pass
