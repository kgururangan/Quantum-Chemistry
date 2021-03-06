function [t1, t2] = build_t1_t2(t1,t2,VM,FM,occ,vir)

    addpath(genpath('/Users/harellab/Dropbox/Hartree Fock/hartree_fock/v4/utils'));

    Nocc = length(occ); Nunocc = length(vir);

    Zocc = ones(Nocc)-eye(Nocc);
    Zunocc = ones(Nunocc)-eye(Nunocc);
    
    foo = FM(occ,occ);
    fov = FM(occ,vir);
    fvv = FM(vir,vir);
    fvo = FM(vir,occ);
    
    Vooov = VM(occ,occ,occ,vir);
    Voovv = VM(occ,occ,vir,vir);
    Vovov = VM(occ,vir,occ,vir);
    Vvovv = VM(vir,occ,vir,vir);
    Vvovo = VM(vir,occ,vir,occ);
    Voovo = VM(occ,occ,vir,occ);
    Vovoo = VM(occ,vir,occ,occ);
    Voooo = VM(occ,occ,occ,occ);
    Vvvvv = VM(vir,vir,vir,vir);
    Vvvvo = VM(vir,vir,vir,occ);
    Vvvoo = VM(vir,vir,occ,occ);

    % Intermediates
    
    chi_mi = foo.*Zocc + ...
             einsum_kg(Vooov,t1,'mnif,fn->mi') + ...
             0.5*einsum_kg(Voovv,t2,'mnef,efin->mi');
    
    chi_me = fov + einsum_kg(Voovv,t1,'mnef,fn->me');
    
    chi_ae = fvv.*Zunocc + ...
             einsum_kg(Vvovv,t1,'anef,fn->ae') - ...
             0.5*einsum_kg(Voovv,t2,'mnef,afmn->ae') - ...
             einsum_kg(chi_me,t1,'me,am->ae');
    
    chit_mi = chi_mi + einsum_kg(chi_me,t1,'me,ei->mi');
    
    chit_anef = Vvovv - 0.5*einsum_kg(Voovv,t1,'mnef,am->anef');
    
    chi_anej = Vvovo - ... 
               0.5*einsum_kg(Voovo,t1,'mnej,am->anej') + ...
               0.5*einsum_kg(chit_anef,t1,'anef,fj->anej');
    
    chi_mnif = Vooov + 0.5*einsum_kg(Voovv,t1,'mnef,ei->mnif');
    
    chi_mbij = Vovoo - ...
                0.5*einsum_kg(Voooo,t1,'mnij,bn->mbij') - ...
                0.5*einsum_kg(chit_anef,t2,'bmef,efij->mbij');
    
    chi_abej = 0.5*Vvvvo - einsum_kg(chi_anej,t1,'anej,bn->abej');
    
    chi_mnij = 0.5*Voooo + ...
               0.25*einsum_kg(Voovv,t2,'mnef,efij->mnij') + ...
               einsum_kg(chi_mnif,t1,'mnif,fj->mnij');
    
    chi_anef = chit_anef - 0.5*einsum_kg(Voovv,t1,'mnef,am->anef');
    
    chit_anej = Vvovo - ... 
                einsum_kg(Voovo,t1,'mnej,am->anej') - ...
                0.5*einsum_kg(Voovv,t2,'mnef,afmj->anej') + ...
                einsum_kg(chi_anef,t1,'anef,fj->anej');
    
    chi_efij = einsum_kg(t1,t1,'ei,fj->efij') + 0.5*t2;
    
    % build T1
    
    TEMP1 = einsum_kg(chi_me,t2,'me,aeim->ai') - ...
            einsum_kg(Vovov,t1,'maie,em->ai') - ...
            einsum_kg(chi_mi,t1,'mi,am->ai') + ...
            einsum_kg(chi_ae,t1,'ae,ei->ai') - ...
            0.5*einsum_kg(Vooov,t2,'mnif,afmn->ai') + ...
            0.5*einsum_kg(Vvovv,t2,'anef,efin->ai');
        
    % build T2

    % standard order
    TEMP2 = einsum_kg(chi_abej,t1,'abej,ei->abij') - ...
            einsum_kg(chit_anej,t2,'anej,ebin->abij') - ...
            0.5*einsum_kg(chi_mbij,t1,'mbij,am->abij') + ...
            0.5*einsum_kg(chi_ae,t2,'ae,ebij->abij') - ...
            0.5*einsum_kg(chit_mi,t2,'mi,abmj->abij') + ...
            0.25*einsum_kg(chi_mnij,t2,'mnij,abmn->abij') + ...
            0.25*einsum_kg(Vvvvv,chi_efij,'abef,efij->abij');
    % (ij)    
    TEMP2_ij = einsum_kg(chi_abej,t1,'abei,ej->abij') - ...
            einsum_kg(chit_anej,t2,'anei,ebjn->abij') - ...
            0.5*einsum_kg(chi_mbij,t1,'mbji,am->abij') + ...
            0.5*einsum_kg(chi_ae,t2,'ae,ebji->abij') - ...
            0.5*einsum_kg(chit_mi,t2,'mj,abmi->abij') + ...
            0.25*einsum_kg(chi_mnij,t2,'mnji,abmn->abij') + ...
            0.25*einsum_kg(Vvvvv,chi_efij,'abef,efji->abij');
    % (ab)    
    TEMP2_ab = einsum_kg(chi_abej,t1,'baej,ei->abij') - ...
            einsum_kg(chit_anej,t2,'bnej,eain->abij') - ...
            0.5*einsum_kg(chi_mbij,t1,'maij,bm->abij') + ...
            0.5*einsum_kg(chi_ae,t2,'be,eaij->abij') - ...
            0.5*einsum_kg(chit_mi,t2,'mi,bamj->abij') + ...
            0.25*einsum_kg(chi_mnij,t2,'mnij,bamn->abij') + ...
            0.25*einsum_kg(Vvvvv,chi_efij,'baef,efij->abij');
    
    % (ab)(ij)    
    TEMP2_abij = einsum_kg(chi_abej,t1,'baei,ej->abij') - ...
            einsum_kg(chit_anej,t2,'bnei,eajn->abij') - ...
            0.5*einsum_kg(chi_mbij,t1,'maji,bm->abij') + ...
            0.5*einsum_kg(chi_ae,t2,'be,eaji->abij') - ...
            0.5*einsum_kg(chit_mi,t2,'mj,bami->abij') + ...
            0.25*einsum_kg(chi_mnij,t2,'mnji,bamn->abij') + ...
            0.25*einsum_kg(Vvvvv,chi_efij,'baef,efji->abij');
        
        
        
    t1 = fvo + TEMP1;
    t2 = Vvvoo + TEMP2 - TEMP2_ij - TEMP2_ab + TEMP2_abij;
            

            
end