import numpy as np

def orthomat(S,tol,kind):
    evalS, U = np.linalg.eigh(S)
    diagS_minushalf = np.diag(evalS**(-0.5))
    if kind == 'symmetric':
        X0 = np.dot(U,np.dot(diagS_minushalf,U.T))
    else:
        X0 = np.dot(U,diagS_minushalf)

    idX = [i for i in range(S.shape[0])]
    for i in range(S.shape[0]):
        for j in range(i+1,S.shape[0]):
            if np.abs(evalS[i]-evalS[j]) < tol:
                J = idX.index(j)
                idX.remove(J)
    X = X0[:,idX]

    return X


def diis_pulay_solver(X_list,diis_resid_list):

    B_dim = len(X_list) + 1
    B = np.empty((B_dim, B_dim))
    B[-1, :] = -1
    B[:, -1] = -1
    B[-1, -1] = 0
    for i in range(len(X_list)):
        for j in range(i,len(X_list)):
            B[i, j] = np.einsum('ij,ij->', diis_resid_list[i], diis_resid_list[j], optimize=True)
            B[j, i] = B[i, j]

    # Build RHS of Pulay equation
    rhs = np.zeros((B_dim))
    rhs[-1] = -1

    # Solve Pulay equation for c_i's
    if np.linalg.det(B) != 0:
        coeff = np.linalg.solve(B, rhs)
        # Build DIIS Fock matrix
        X = np.zeros_like(X_list[0])
        for x in range(B_dim - 1):
            X += coeff[x] * X_list[x]
        return X
    else:
        return X_list[-1]
