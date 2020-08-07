function [HBar] = build_HBar( t1, t2, sys )

    addpath(genpath('/Users/harellab/Dropbox/Hartree Fock/hartree_fock/v4/utils'));

    tic
    
    Nunocc = sys.Nunocc; Nocc = sys.Nocc;
    
    tau = t2 + 0.5*(einsum_kg(t1,t1,'ai,bj->abij') - ...
                    einsum_kg(t1,t1,'aj,bi->abij') - ...
                    einsum_kg(t1,t1,'bi,aj->abij') + ...
                    einsum_kg(t1,t1,'bj,ai->abij')); 
    
    taubar = t2 + 0.25*(einsum_kg(t1,t1,'ai,bj->abij') - ...
                        einsum_kg(t1,t1,'aj,bi->abij') - ...
                        einsum_kg(t1,t1,'bi,aj->abij') + ...
                        einsum_kg(t1,t1,'bj,ai->abij')); 

    Fbar_ae =   sys.fvv - ...
                0.5*einsum_kg(sys.fov,t1,'me,am->ae') + ...
                einsum_kg(sys.Vvovv,t1,'amef,fm->ae') - ...
                0.5*einsum_kg(taubar,sys.Voovv,'afmn,mnef->ae');
                
    Fbar_mi = sys.foo + ...
              0.5*einsum_kg(sys.fov,t1,'me,ei->mi')  + ...
              einsum_kg(sys.Vooov,t1,'mnie,en->mi') + ...
              0.5*einsum_kg(sys.Voovv,taubar,'mnef,efin->mi');
                
    Fbar_me = sys.fov + einsum_kg(sys.Voovv,t1,'mnef,fn->me');
    
    Wbar_mnij = sys.Voooo + ...
                einsum_kg(sys.Vooov,t1,'mnie,ej->mnij') - einsum_kg(sys.Vooov,t1,'mnje,ei->mnij') + ...
                0.25*einsum_kg(sys.Voovv,tau,'mnef,efij->mnij');
        
    Wbar_abef = sys.Vvvvv - ...
                einsum_kg(sys.Vvovv,t1,'amef,bm->abef') + einsum_kg(sys.Vvovv,t1,'bmef,am->abef') + ...
                0.25*einsum_kg(sys.Voovv,tau,'mnef,abmn->abef');
            
    Wbar_mbej = sys.Vovvo + ...
                einsum_kg(sys.Vovvv,t1,'mbef,fj->mbej') - ...
                einsum_kg(sys.Voovo,t1,'mnej,bn->mbej') - ...
                0.5*einsum_kg(sys.Voovv,t2,'mnef,fbjn->mbej') - ...
                einsum_kg(einsum_kg(sys.Voovv,t1,'mnef,fj->mnej'),t1,'mnej,bn->mbej');
            
    % Hbar 1-body components
    H1 = cell(2,2);
    
    H1{2,2} = Fbar_ae - 0.5*einsum_kg(Fbar_me,t1,'me,am->ae');
    
    H1{1,1} = Fbar_mi + 0.5*einsum_kg(Fbar_me,t1,'me,ei->mi');
    
    H1{1,2} = Fbar_me;
    
    H1{2,1} = zeros(Nunocc,Nocc); % zero by CCSD
    
    % Hbar 2-body components
    H2 = cell(2,2,2,2);
    
    H2{1,1,2,2} = sys.Voovv;
    
    H2{1,1,1,1} = Wbar_mnij + 0.25*einsum_kg(sys.Voovv,tau,'mnef,efij->mnij'); % W_mnij
    
    H2{2,2,2,2} = Wbar_abef + 0.25*einsum_kg(sys.Voovv,tau,'mnef,abmn->abef'); % W_abef
    
    H2{1,2,2,1} = Wbar_mbej - 0.5*einsum_kg(sys.Voovv,t2,'mnef,fbjn->mbej'); % W_mbej
    
    H2{2,1,1,2} = permute(H2{1,2,2,1},[2,1,4,3]);
    
    H2{2,1,2,1} = -permute(H2{1,2,2,1},[2,1,3,4]);
    
    H2{1,1,1,2} = sys.Vooov + einsum_kg(sys.Voovv,t1,'mnfe,fi->mnie'); % W_mnie
    
    H2{1,1,2,1} = -permute(H2{1,1,1,2},[1,2,4,3]);
    
    H2{2,1,2,2} = sys.Vvovv - einsum_kg(sys.Voovv,t1,'nmef,an->amef'); % W_amef
    
    H2{1,2,2,2} = -permute(H2{2,1,2,2},[2,1,3,4]);
    
    H2{1,2,1,1} =   sys.Vovoo - einsum_kg(H1{1,2},t2,'me,beij->mbij') - einsum_kg(H2{1,1,1,1},t1,'mnij,bn->mbij') + ...
                    0.5*einsum_kg(sys.Vovvv,tau,'mbef,efij->mbij') + ...
                    einsum_kg(sys.Vooov,t2,'mnie,bejn->mbij') - einsum_kg(sys.Vooov,t2,'mnje,bein->mbij') + ...
                    einsum_kg(sys.Vovvo,t1,'mbej,ei->mbij') - einsum_kg(sys.Vovvo,t1,'mbei,ej->mbij') - ...
                    einsum_kg(einsum_kg(sys.Voovv,t2,'mnef,bfnj->mbej'),t1,'mbej,ei->mbij') + ...
                    einsum_kg(einsum_kg(sys.Voovv,t2,'mnef,bfni->mbei'),t1,'mbei,ej->mbij'); % W_mbij

    H2{2,1,1,1} = -permute(H2{1,2,1,1},[2,1,3,4]);
                
    H2{2,2,2,1} =   sys.Vvvvo - einsum_kg(H1{1,2},t2,'me,abmi->abei') + einsum_kg(H2{2,2,2,2},t1,'abef,fi->abei') + ...
                    0.5*einsum_kg(sys.Voovo,tau,'mnei,abmn->abei') - ...
                    einsum_kg(sys.Vovvv,t2,'mbef,afmi->abei') + einsum_kg(sys.Vovvv,t2,'maef,bfmi->abei') - ...
                    einsum_kg(sys.Vovvo,t1,'mbei,am->abei') + einsum_kg(sys.Vovvo,t1,'maei,bm->abei') + ...
                    einsum_kg(einsum_kg(sys.Voovv,t2,'mnef,bfni->mbei'),t1,'mbei,am->abei') - ...
                    einsum_kg(einsum_kg(sys.Voovv,t2,'mnef,afni->maei'),t1,'maei,bm->abei'); % W_abei

    H2{2,2,1,2} = -permute(H2{2,2,2,1},[1,2,4,3]);
                
    H2{2,2,1,1} = zeros(Nunocc,Nunocc,Nocc,Nocc); % zero by CCSD
                
    % Hbar 3-body matrix elements
    H3 = cell(2,2,2,2,2,2);
     
    H3{2,2,1,1,1,2} =  -einsum_kg(sys.Voovo,t2,'mnej,abin->abmije') +...
                        einsum_kg(sys.Voovo,t2,'mnei,abjn->abmije') - ...
                        einsum_kg(einsum_kg(sys.Voovv,t1,'mnef,fi->mnei'),t2,'mnei,abnj->abmije') + ...
                        einsum_kg(einsum_kg(sys.Voovv,t1,'mnef,fj->mnej'),t2,'mnej,abni->abmije') + ...
                        einsum_kg(sys.Vovvv,t2,'mbef,afij->abmije') - ...
                        einsum_kg(sys.Vovvv,t2,'maef,bfij->abmije') - ...
                        einsum_kg(einsum_kg(sys.Voovv,t1,'mnef,an->maef'),t2,'maef,fbij->abmije') + ...
                        einsum_kg(einsum_kg(sys.Voovv,t1,'mnef,bn->mbef'),t2,'mbef,faij->abmije'); % W_abmije

    H3{2,1,1,1,1,2} = einsum_kg(sys.Voovv,t2,'mnfe,afij->amnije'); % W_amnije
    
    H3{2,2,1,1,2,2} = -einsum_kg(sys.Voovv,t2,'mnfe,abin->abmief'); % W_abmief
    
    H3{2,1,1,1,2,1} = -permute(H3{2,1,1,1,1,2},[1,2,3,4,6,5]);
    
    H3{2,1,2,1,2,2} = -permute(H3{2,2,1,1,2,2},[1,3,2,4,5,6]);
    
    HBar = {H1, H2, H3};
    
    
    
    fprintf('\nHBar succesfully built in %4.2f seconds\n',toc);

end

