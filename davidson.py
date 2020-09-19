import numpy as np
import argparse

def gramschmidt(X):

    d, n = X.shape
    m = min(d,n)
    R = np.zeros((m,n))
    Q = np.zeros((d,m))
    Q2 = np.zeros((d,n))

    for i in range(m):
        v = X[:,i]
        for j in range(i):
            R[j,i] = np.dot(Q[:,j].T,v)
            v -= R[j,i]*Q[:,j]
        R[i,i] = np.linalg.norm(v)
        Q[:,i] = v/R[i,i]
    R[:,m+1:] = np.dot(Q.T,X[:,m+1:])

    return Q,R

def diagonal_guess(nroot,D):
    idx = np.argsort(D)
    I = np.eye(len(D))
    return I[:,idx[:nroot]]

def update_R(r,e,D):
    for i in range(len(r)):
        r[i] *= 1/(e - D[i])
    return r

def orthogonalize_root(r,B):
    for i in range(B.shape[1]):
        b = B[:,i]/np.linalg.norm(B[:,i])
        r -= np.dot(b.T,r)*b
    return r

def davidson(A,D,nroot,nvec,maxit,tol,thresh_vec):

    vec_dim = A.shape[0]
    max_size = nroot*nvec

    B = diagonal_guess(nroot,D)
    curr_size = nroot

    for iter in range(maxit):

        B,_ = gramschmidt(B)

        sigma = np.dot(A,B)

        G = np.dot(B.T,sigma)
        
        # MUST USE np.linalg.eig NOT np.linalg.eigh
        eval_E, alpha = np.linalg.eig(G)

        idx = np.argsort(eval_E)
        eigval = eval_E[idx[:nroot]]
        alpha = alpha[:,idx[:nroot]]
        V = np.dot(B,alpha)

        add_B = np.zeros((vec_dim,nroot))
        resid_norm = np.zeros(nroot)
        ct_add = 0
        print('\nIter - {}    Subspace Dim = {}'.format(iter+1,curr_size))
        print('----------------------------------------------------')
        for j in range(nroot):
            r = np.dot(sigma,alpha[:,j]) - eigval[j]*V[:,j]
            resid_norm[j] = np.linalg.norm(r)
            q = update_R(r,eigval[j],D)
            q = q/np.linalg.norm(q)
            if ct_add > 0:
                q_orth = orthogonalize_root(q,np.concatenate((B,add_B[:,:ct_add]),axis=1))
            else:
                q_orth = orthogonalize_root(q,B)
            if np.linalg.norm(q_orth) > thresh_vec:
                add_B[:,ct_add] = q_orth/np.linalg.norm(q_orth)
                ct_add += 1

            print('Root - {}      e = {:0.8f}      |r| = {:0.8f}'.format(j+1,eigval[j],resid_norm[j]))

        if all(resid_norm <= tol):
            print('\nDavidson successfully converged!')
            break
        else:
            if curr_size >= max_size:
                print('Restarting and collapsing')
                B = np.dot(B,alpha)
                curr_size = nroot
            else:
                B = np.concatenate((B,np.asarray(add_B[:,:ct_add])),axis=1)
                curr_size = B.shape[1]

    return eigval, V
    
def main(args):
    
    n = args.size
    nvec = args.nvec
    tol = args.tolerance
    nroot = args.roots
    maxit = args.maxit
    thresh_vec = args.threshold

    sparsity = 10
    hermiticity = 0.5

    R = np.random.rand(n,n)
    A = np.diag(np.arange(1,n+1)) + sparsity*( hermiticity*R + (1.0-hermiticity)*R.T )
    D = np.diagonal(A)


    E_dav = np.zeros(nroot)
    V_dav = np.zeros((n,nroot))
    E_dav, V_dav = davidson(A,D,nroot,nvec,maxit,tol,thresh_vec)


    E_exact, V_exact = np.linalg.eig(A)
    idx = np.argsort(E_exact)
    E_exact = [E_exact[i] for i in idx]
    print('\n Root      Davidson         Exact')
    print('====================================')
    for i in range(nroot):
        print('  {}      {:0.8f}      {:0.8f}'.format(i+1,E_dav[i],E_exact[i]))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-v','--nvec',type=int,help='Max number of vectors per root',default=5)
    parser.add_argument('-t','--tolerance',type=float,help='Eigenvector residuum tolerance',default=1e-07)
    parser.add_argument('-n','--size',type=int,help='N x N size of test matrix',default=100)
    parser.add_argument('-r','--roots',type=int,help='Number of roots',default=5)
    parser.add_argument('-m','--maxit',type=int,help='Maximum number of iterations',default=30)
    parser.add_argument('-th','--threshold',type=float,help='Subspace addition threshold',default=1e-04)
    parser.add_argument('-s','--sparsity',type=float,help='Sparsity of random test matrix',default=2)
    parser.add_argument('-y','--hermitian',type=float,help='Degree of hermiticity of test matrix',default=0.5)
    args = parser.parse_args()
    main(args)

