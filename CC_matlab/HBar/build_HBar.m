function [HBar] = build_HBar( t1, t2, sys )

    addpath(genpath('/Users/harellab/Dropbox/Hartree Fock/hartree_fock/v4/utils'));
    
    VM = sys.VM; FM = sys.FM;

    tic
    
    Nocc = size(t1,2); Nunocc = size(t1,1); Norb = size(FM,1);
    iocc = 1:Nocc; ivir = Nocc+1:Norb;

%    fvo = FM(ivir, iocc);
    fov = FM(iocc, ivir);
    foo = FM(iocc, iocc);
    fvv = FM(ivir, ivir);
    
%    Vovvv = VM(iocc,ivir,ivir,ivir);
    Voovv = VM(iocc,iocc,ivir,ivir);
    Vooov = VM(iocc,iocc,iocc,ivir);
    Vovvv = VM(iocc,ivir,ivir,ivir);
%    Vovov = VM(iocc,ivir,iocc,ivir);
%    Vvoov = VM(ivir,iocc,iocc,ivir);
%    Vvovo = VM(ivir,iocc,ivir,iocc);
    Voovo = VM(iocc,iocc,ivir,iocc);
    Vvovv = VM(ivir,iocc,ivir,ivir);
    Vvvvv = VM(ivir,ivir,ivir,ivir);
%    Vvooo = VM(ivir,iocc,iocc,iocc);
    Voooo = VM(iocc,iocc,iocc,iocc);
%    Vvvov = VM(ivir,ivir,iocc,ivir);
    Vvvvo = VM(ivir,ivir,ivir,iocc);
    Vovoo = VM(iocc,ivir,iocc,iocc);
    Vovvo = VM(iocc,ivir,ivir,iocc);
    
    tau = t2 + 0.5*(einsum(t1,t1,'ai,bj->abij') - ...
                    einsum(t1,t1,'aj,bi->abij') - ...
                    einsum(t1,t1,'bi,aj->abij') + ...
                    einsum(t1,t1,'bj,ai->abij')); 
    
    taubar = t2 + 0.25*(einsum(t1,t1,'ai,bj->abij') - ...
                        einsum(t1,t1,'aj,bi->abij') - ...
                        einsum(t1,t1,'bi,aj->abij') + ...
                        einsum(t1,t1,'bj,ai->abij')); 

    Fbar_ae =   fvv - ...
                0.5*einsum(fov,t1,'me,am->ae') + ...
                einsum(Vvovv,t1,'amef,fm->ae') - ...
                0.5*einsum(taubar,Voovv,'afmn,mnef->ae');
                
    Fbar_mi = foo + ...
              0.5*einsum(fov,t1,'me,ei->mi')  + ...
              einsum(Vooov,t1,'mnie,en->mi') + ...
              0.5*einsum(Voovv,taubar,'mnef,efin->mi');
                
    Fbar_me = fov + einsum(Voovv,t1,'mnef,fn->me');
    
    Wbar_mnij = Voooo + ...
                einsum(Vooov,t1,'mnie,ej->mnij') - einsum(Vooov,t1,'mnje,ei->mnij') + ...
                0.25*einsum(Voovv,tau,'mnef,efij->mnij');
        
    Wbar_abef = Vvvvv - ...
                einsum(Vvovv,t1,'amef,bm->abef') + einsum(Vvovv,t1,'bmef,am->abef') + ...
                0.25*einsum(Voovv,tau,'mnef,abmn->abef');
            
    Wbar_mbej = Vovvo + ...
                einsum(Vovvv,t1,'mbef,fj->mbej') - ...
                einsum(Voovo,t1,'mnej,bn->mbej') - ...
                0.5*einsum(Voovv,t2,'mnef,fbjn->mbej') - ...
                einsum(einsum(Voovv,t1,'mnef,fj->mnej'),t1,'mnej,bn->mbej');
            
    % Hbar 1-body components
    H1 = cell(2,2);
    
    H1{2,2} = Fbar_ae - 0.5*einsum(Fbar_me,t1,'me,am->ae');
    
    H1{1,1} = Fbar_mi + 0.5*einsum(Fbar_me,t1,'me,ei->mi');
    
    H1{1,2} = Fbar_me;
    
    H1{2,1} = zeros(Nunocc,Nocc); % zero by CCSD
    
    % Hbar 2-body components
    H2 = cell(2,2,2,2);
    
    H2{1,1,2,2} = Voovv;
    
    H2{1,1,1,1} = Wbar_mnij + 0.25*einsum(Voovv,tau,'mnef,efij->mnij'); % W_mnij
    
    H2{2,2,2,2} = Wbar_abef + 0.25*einsum(Voovv,tau,'mnef,abmn->abef'); % W_abef
    
    H2{1,2,2,1} = Wbar_mbej - 0.5*einsum(Voovv,t2,'mnef,fbjn->mbej'); % W_mbej
    
    H2{1,1,1,2} = Vooov + einsum(Voovv,t1,'mnfe,fi->mnie'); % W_mnie
    
    H2{2,1,2,2} = Vvovv - einsum(Voovv,t1,'nmef,an->amef'); % W_amef
    
    H2{1,2,1,1} =   Vovoo - einsum(H1{1,2},t2,'me,beij->mbij') - einsum(H2{1,1,1,1},t1,'mnij,bn->mbij') + ...
                    0.5*einsum(Vovvv,tau,'mbef,efij->mbij') + ...
                    einsum(Vooov,t2,'mnie,bejn->mbij') - einsum(Vooov,t2,'mnje,bein->mbij') + ...
                    einsum(Vovvo,t1,'mbej,ei->mbij') - einsum(Vovvo,t1,'mbei,ej->mbij') - ...
                    einsum(einsum(Voovv,t2,'mnef,bfnj->mbej'),t1,'mbej,ei->mbij') + ...
                    einsum(einsum(Voovv,t2,'mnef,bfni->mbei'),t1,'mbei,ej->mbij'); % W_mbij
                
    H2{2,2,2,1} =   Vvvvo - einsum(H1{1,2},t2,'me,abmi->abei') + einsum(H2{2,2,2,2},t1,'abef,fi->abei') + ...
                    0.5*einsum(Voovo,tau,'mnei,abmn->abei') - ...
                    einsum(Vovvv,t2,'mbef,afmi->abei') + einsum(Vovvv,t2,'maef,bfmi->abei') - ...
                    einsum(Vovvo,t1,'mbei,am->abei') + einsum(Vovvo,t1,'maei,bm->abei') + ...
                    einsum(einsum(Voovv,t2,'mnef,bfni->mbei'),t1,'mbei,am->abei') - ...
                    einsum(einsum(Voovv,t2,'mnef,afni->maei'),t1,'maei,bm->abei'); % W_abei
                
    H2{2,2,1,1} = zeros(Nunocc,Nunocc,Nocc,Nocc); % zero by CCSD
                
    % Hbar 3-body matrix elements
    H3 = cell(2,2,2,2,2,2);
     
    H3{2,2,1,1,1,2} =  -einsum(Voovo,t2,'mnej,abin->abmije') +...
                        einsum(Voovo,t2,'mnei,abjn->abmije') - ...
                        einsum(einsum(Voovv,t1,'mnef,fi->mnei'),t2,'mnei,abnj->abmije') + ...
                        einsum(einsum(Voovv,t1,'mnef,fj->mnej'),t2,'mnej,abni->abmije') + ...
                        einsum(Vovvv,t2,'mbef,afij->abmije') - ...
                        einsum(Vovvv,t2,'maef,bfij->abmije') - ...
                        einsum(einsum(Voovv,t1,'mnef,an->maef'),t2,'maef,fbij->abmije') + ...
                        einsum(einsum(Voovv,t1,'mnef,bn->mbef'),t2,'mbef,faij->abmije'); % W_abmije

    H3{2,1,1,1,1,2} = einsum(Voovv,t2,'mnfe,afij->amnije'); % W_amnije
    
    H3{2,2,1,1,2,2} = -einsum(Voovv,t2,'mnfe,abin->abmief'); % W_abmief
    
    HBar = {H1, H2, H3};
    
    
    
    fprintf('\nHBar succesfully built in %4.2f seconds\n',toc);

end

