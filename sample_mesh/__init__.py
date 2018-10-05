from .py_sample_mesh import *

def read_obj(filename):
    V = []
    F = []
    N = [[]]
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('v'):
                V.append([float(vi) for vi in line.split()[1:]])
            elif line.startswith('f'):
                F.append([int(fi.split('/')[0])-1 for fi in line.split()[1:]])

    V = np.array(V)
    F = np.array(F)

    return V, F
