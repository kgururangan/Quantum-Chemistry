function [Ecc] = cc_energy(t1,t2,VM,FM,occ,unocc)

    %addpath(genpath('/Users/harellab/Dropbox/Hartree Fock/hartree_fock/v4/utils'));
    
    fov = FM(occ,unocc);
    Voovv = VM(occ,occ,unocc,unocc);
    
    Ecc = einsum(fov,t1,'ia,ai->') + ...
          0.25*einsum(Voovv,t2,'ijab,abij->') + ...
          0.5*einsum(einsum(Voovv,t1,'ijab,ai->jb'),t1,'jb,bj->');

end