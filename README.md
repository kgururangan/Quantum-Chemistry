# Quantum-Chemistry
Collection of codes used for quantum chemistry calculations based on Hartree-Fock and Coupled Cluster (CC) methodologies. These are all my own codes developed for learning and testing purposes.

Contents of Repository:

CC_matlab - most up-to-date set of coupled-cluster codes. Available calculations:

    - CCSD (spinorbital & spin-integrated)
    - CR-CC(2,3) (spinorbital)
    - CCSDT (spinorbital)
    - Active-space CCSDT(III) (spinorbital)
    - Left-CCSD (spinorbital)
    - EOM-CCSD (spinorbital)
    - Left-EOMCCSD (spinorbital)
  
  
Routines can be run using one- and two-body molecular orbital integrals inputs (along with system parameters like number of occupied and unoccupied alpha/beta electrons) which can be found in CC_matlab_tests. Run routines as shown in main.m

Numerov - Matlab implementation of the Numerov-Cooley algorithm for numerically integrating the 1D Schrodinger equation. The main.m script applies this method to solve for the 1D vibrational eigenstates of the H2 stretching motion using a simple restricted Hartree-Fock (RHF) potential energy surface.

CT_Hamiltonian - collection of codes used for crude modeling of charge-transfer processes in dimers

Davidson_matlab_test - Matlab implementation of davidson diagonalization routines. Used in EOMCC matlab codes

gaussian_integrals - Fortran implementation of McMurchie-Davidson atomic orbital integral evaluation scheme. I largely translated the Python implementation originally written by Joshua Goings https://github.com/jjgoings/McMurchie-Davidson

hf_python_v* - My initial Python implementation of Hartree-Fock (with integral routines directly taken from Joshua Goings again) which lead into spinorbital CCSD and EOMCCSD. Big and unwieldy code which is streamlined and properly implemented in CC_matlab (however, it's useful for generated one- and two-body integrals needed to run the Matlab version).

Davidson_test - Python implementation of Davidson routine (use Matlab versions as they are more stable)

Permutation Sign - Useful function for evaluating the sign of a permutation

Tensor Contraction - An exercise in unravelling the amazing Numpy einsum function and writing my own einsum from scratch

Wicks_Theorem_v2 - An attempt at creating a Wick's Theorem evaluation engine for automated equation derivation in coupled-cluster methods. Not exactly complete (or fast) but can evaluate fully-contracted second-quantized expressions and classify the results in Hugenholtz diagrams

