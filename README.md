# Quantum-Chemistry
Collection of codes used for quantum chemistry calculations based on Hartree-Fock and Coupled Cluster (CC) methodologies. These are all my own codes developed for learning and testing purposes.

Contents of Repository:

**CC_matlab**

Spin-integrated coupled-cluster codes written in Matlab that implement coupled-cluster (CC), completely renormalized (CR) CC, equation-of-motion (EOM) CC, and active-space CC methods. All routines are vectorized and equations are highly factorized using intermediates derived from the similarity-transformed Hamiltonian. The code is meant to serve as a pedagogical template for coding coupled-cluster methods as well as a sandbox for testing out new ideas in a simple coding language.

The codes are simply run using one-body and two-body molecular orbital integrals obtained from, for instance, a quantum chemistry software such as GAMESS or Quantum Package 2.0. It is compatible with RHF, UHF, and ROHF reference states. Also contains SCF solvers so that molecular orbital integrals can be generated in a self-contained manner, although there is currently no symmetry adaptation. Currently supported calculations include:

    - RHF + RHF Analytical Gradients
    - CCSD + Left-CCSD
    - EOMCCSD + Left-EOMCCSD
    - CR-CC(2,3) + CR-EOMCC(2,3)
    - CCSDT
    - Active-space CCSDt variants II and III
    - Moment-corrected CC(t;3) II and III
  
Routines can be run using one- and two-body molecular orbital integrals inputs (along with system parameters like number of occupied and unoccupied alpha/beta electrons) which can be found in CC_matlab_tests. Run routines as shown in main.m

**CCpy**

A Python software suite of coupled-cluster codes, similar to the contents of CC_matlab, including CC, CR-CC, and EOMCC methods. It is based on the optimized Numpy einsum routine as well as extremeley fast updates written in Fortran and called within Python using ```f2py```. Includes the methods in CC_matlab as well as CR-CC(2,4).

**einsum_f90**

A first attempt at writing an Einstein summation function in Modern Fortran.

**Numerov**

Matlab implementation of the Numerov-Cooley algorithm for numerically integrating the 1D Schrodinger equation. The main.m script applies this method to solve for the 1D vibrational eigenstates of the H2 stretching motion using a simple restricted Hartree-Fock (RHF) potential energy surface.

**CT_Hamiltonian**

Collection of codes used for modeling of charge-transfer processes in dimers.

**Davidson_matlab_test**

Matlab implementation of davidson diagonalization routines. Used in EOMCC matlab codes

**gaussian_integrals**

Fortran implementation of McMurchie-Davidson atomic orbital integral evaluation scheme. I largely translated the Python implementation originally written by Joshua Goings https://github.com/jjgoings/McMurchie-Davidson

**hf_python_v **

My initial Python implementation of Hartree-Fock (with integral routines directly taken from Joshua Goings again) which lead into spinorbital CCSD and EOMCCSD. Big and unwieldy code which is streamlined and properly implemented in CC_matlab.

**Davidson_test**

Python implementation of Davidson routine (use Matlab versions as they are more stable)

**Permutation Sign**

Useful function for evaluating the sign of a permutation

**Tensor Contraction**

An exercise in unravelling the amazing Numpy einsum function and writing my own einsum from scratch

**Wicks_Theorem_v2**

An attempt at creating a Wick's Theorem evaluation engine for automated equation derivation in coupled-cluster methods. Not exactly complete (or fast) but can evaluate fully-contracted second-quantized expressions and classify the results in Hugenholtz diagrams

