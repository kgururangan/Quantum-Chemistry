import numpy as np
import argparse 
from matplotlib import pyplot as plt


def get_pars(E,V0,m):

    k_sq = 2*m*np.abs(E)
    kappa_sq = 2*m*np.abs(E-V0)

    if E < 0:
        k = 1j*np.sqrt(k_sq)
    else:
        k = np.sqrt(k_sq)

    if E > V0:
        kappa = np.sqrt(kappa_sq)
    else:
        kappa = 1j*np.sqrt(kappa_sq)

    return k, kappa

def get_system(k,kappa,a):

    A = np.array([[1.0, -1.0, -1.0, 0.0],
                  [-k, -kappa, kappa, 0.0],
                  [0.0, np.exp(1j*kappa*a), np.exp(-1j*kappa*a), -np.exp(1j*k*a)],
                  [0.0, kappa*np.exp(1j*kappa*a), -kappa*np.exp(-1j*kappa*a), -k*np.exp(1j*k*a)]], dtype=np.complex64)

    b = np.array([-1.0, -k, 0.0, 0.0],dtype=np.complex64)
    

    return A, b

def solve_system(A,b):

    x = np.linalg.solve(A,b)
    R = np.real(np.sqrt(np.conj(x[0])*x[0]))
    T = np.real(np.sqrt(np.conj(x[3])*x[3]))

    return R, T

def main():

    V0 = 1
    a = 1.0
    m = 20
    E = np.linspace(1.1,10,1000)

    T = np.zeros(len(E))
    R = np.zeros(len(E))

    for i in range(len(E)):

        print('E = {}'.format(E[i]))

        k, kappa = get_pars(E[i], V0, m) 
        A, b = get_system(k, kappa, a)
        R[i], T[i] = solve_system(A, b)

        #print('R = {}'.format(r**2))
        #print('T = {}'.format(t**2))
                

    plt.plot(E/V0,T*100.0,label='Transmission')
    plt.plot(E/V0,R*100.0,label='Reflection')
    plt.xlabel(r'$E / V_0$')
    plt.ylabel('[%]')
    plt.legend(loc=5)
    plt.show()
    plt.close()


if __name__ == '__main__':

    main()
