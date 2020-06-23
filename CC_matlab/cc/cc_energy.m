function [Ecc] = cc_energy(t1,t2,sys)

    %addpath(genpath('/Users/harellab/Dropbox/Hartree Fock/hartree_fock/v4/utils'));

    
    Ecc = einsum_kg(sys.fov,t1,'ia,ai->') + ...
          0.25*einsum_kg(sys.Voovv,t2,'ijab,abij->') + ...
          0.5*einsum_kg(einsum_kg(sys.Voovv,t1,'ijab,ai->jb'),t1,'jb,bj->');

end