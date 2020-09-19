import numpy as np
import argparse

def diis_xtrap(x_list, diis_resid):

    vec_dim, diis_dim = np.shape(x_list)
    B_dim = diis_dim + 1
    B = -1.0*np.ones((B_dim,B_dim))

    for i in range(diis_dim):
        for j in range(diis_dim):
            B[i,j] = np.dot(diis_resid[:,i].T,diis_resid[:,j])
    B[-1,-1] = 0.0

    rhs = np.zeros(B_dim)
    rhs[-1] = -1.0

    coeff = np.linalg.solve(B,rhs)
    x_xtrap = np.zeros(vec_dim)
    for i in range(diis_dim):
        x_xtrap += coeff[i]*x_list[:,i]

    return x_xtrap

def main(args):

    n = args.matrix_size
    maxiter = args.maxit
    omega = args.omega
    diis_size = args.diis_size
    tol = args.tolerance
    sparsity = args.sparsity
    flag_diis = args.diis

    A = np.zeros((n,n))
    b = np.random.rand(n)

    print('Building matrix...')
    for i in range(n):
        for j in range(n):
            if i == j:
                A[i,i] = float(i+1) + sparsity*np.random.rand(1)
            else:
                A[i,j] = sparsity*np.random.rand(1)

    x = b
    Az = A*(np.ones((n,n)) - np.eye(n))
    x_xtrap = np.zeros(n)
    x_list = np.zeros((n,diis_size))
    diis_resid = np.zeros((n,diis_size))

    res = np.dot(A,x) - b

    print('Jacobi solver...')
    it = 0
    while np.linalg.norm(res) > tol and it < maxiter:

        u = b - np.dot(Az,x)
        for j in range(n):
            x[j] = (1.0-omega)*x[j] + omega*u[j]/A[j,j]

        res = np.dot(A,x) - b

        if flag_diis:
            diis_resid[:,np.mod(it,diis_size)] = res
            x_list[:,np.mod(it,diis_size)] = x
            if it > diis_size+1:
                x = diis_xtrap(x_list,diis_resid)

        print('Iter - {}    Residuum = {}'.format(it,np.linalg.norm(res)))

        it += 1

    #print('|x - x_exact| = {}'.format(np.linalg.norm(x-x_exact)))

    return



if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-n','--matrix_size',type=int,help='Size of matrix',default=2000)
    parser.add_argument('-m','--maxit',type=int,help='Maximum number of jacobi iterations',default=100)
    parser.add_argument('-w','--omega',type=np.float64,help='Weighted jacobi damping parameter',default=1.0)
    parser.add_argument('-d','--diis_size',type=int,help='DIIS subspace dimension',default=5)
    parser.add_argument('-t','--tolerance',type=np.float64,help='Convergence tolerance',default=1e-08)
    parser.add_argument('-s','--sparsity',type=np.float64,help='Matrix sparsity',default=1e-03)
    parser.add_argument('--diis',type=bool,help='Boolean indicating whether DIIS is used',default=False)
    args = parser.parse_args()
    main(args)
