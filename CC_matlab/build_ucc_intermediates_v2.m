function [chi1A, chi1B, chi2A, chi2B, chi2C] = build_ucc_intermediates_v2(t1a, t1b, t2a, t2b, t2c, sys)


% fock masks
Zocc_alpha = ones(sys.Nocc_alpha) - eye(sys.Nocc_alpha);
Zvir_alpha = ones(sys.Nvir_alpha) - eye(sys.Nvir_alpha);
Zocc_beta = ones(sys.Nocc_beta) - eye(sys.Nocc_beta);
Zvir_beta = ones(sys.Nvir_beta) - eye(sys.Nvir_beta);

% 1-body intermediates

chi1A.me = sys.fa_ov + einsum(sys.vA_oovv,t1a,'mnef,fn->me') + einsum(sys.vB_oovv,t1b,'mnef,fn->me');

chi1B.me = sys.fb_ov + einsum(sys.vC_oovv,t1b,'nmfe,fn->me') + einsum(sys.vB_oovv,t1a,'nmfe,fn->me');

chi1A.mi = sys.fa_oo.*Zocc_alpha +...
         einsum(sys.vA_ooov,t1a,'mnif,fn->mi') +...
         einsum(sys.vB_ooov,t1b,'mnif,fn->mi') +...
         0.5*einsum(sys.vA_oovv,t2a,'mnef,efin->mi') + ...
         einsum(sys.vB_oovv,t2b,'mnef,efin->mi');

chi1A.mi_bar = chi1A.mi + 0.5*einsum(sys.fa_ov,t1a,'me,ei->mi');
     
chi1B.mj = sys.fb_oo.*Zocc_beta +...
         einsum(sys.vC_ooov,t1b,'mnje,en->mj') +...
         einsum(sys.vB_oovo,t1a,'nmej,en->mj') +...
         0.5*einsum(sys.vC_oovv,t2c,'nmfe,fenj->mj') + ...
         einsum(sys.vB_oovv,t2b,'nmfe,fenj->mj');
     
chi1B.mj_bar = chi1B.mj + 0.5*einsum(sys.fb_ov,t1b,'me,ej->mj');
     
chi1A.ae = sys.fa_vv.*Zvir_alpha + ...
           einsum(sys.vA_vovv,t1a,'anef,fn->ae') +...
           einsum(sys.vB_vovv,t1b,'anef,fn->ae') -...
           0.5*einsum(sys.vA_oovv,t2a,'mnef,afmn->ae') -...
           einsum(sys.vB_oovv,t2b,'mnef,afmn->ae');
       
chi1A.ae_bar = chi1A.ae - 0.5*einsum(sys.fa_ov,t1a,'me,am->ae');
       
chi1B.be = sys.fb_vv.*Zvir_beta + ...
           einsum(sys.vC_ovvv,t1b,'nbfe,fn->be') +...
           einsum(sys.vB_ovvv,t1a,'nbfe,fn->be') -...
           0.5*einsum(sys.vC_oovv,t2c,'nmfe,fbnm->be') -...
           einsum(sys.vB_oovv,t2b,'nmfe,fbnm->be');
       
chi1B.be_bar = chi1B.be - 0.5*einsum(sys.fb_ov,t1b,'me,bm->be');


% 2-body intermediates

% 1p / (2p,1h)
chi2A.amef_bar = sys.vA_vovv - 0.5*einsum(sys.vA_oovv,t1a,'nmef,an->amef');
chi2B.amef_bar = sys.vB_vovv - einsum(sys.vB_oovv,t1a,'nmef,an->amef');
chi2B.maef_bar = sys.vB_ovvv - einsum(sys.vB_oovv,t1b,'mnef,an->maef');
chi2C.abef_bar = sys.vC_vovv - 0.5*einsum(sys.vC_oovv,t1b,'nmef,an->amef');

% 1h / (1p,2h)
chi2A.mnie_bar = sys.vA_ooov + 0.5*einsum(sys.vA_oovv,t1a,'mnfe,fi->mnie');
chi2B.mnie_bar = sys.vA_ooov + einsum(sys.vB_oovv,t1a,'mnfe,fi->mnie');
chi2B.nmei_bar = sys.vA_oovo + einsum(sys.vB_oovv,t1b,'nmef,fi->nmei');
chi2C.mnie_bar = sys.vA_ooov + 0.5*einsum(sys.vC_oovv,t1a,'mnfe,fi->mnie');

% 2p / 2p
chi2A.abef

% 2h / 2h





% chiC                  
chi2C.bmje = sys.vC_voov +...
             einsum(sys.vC_vovv,t1b,'bmfe,fj->bmje') -...
             einsum(sys.vC_oovo,t1b,'mnej,bn->bmje') +...
             einsum(einsum(sys.vC_oovv,t1b,'mnef,bn->mbef'),t1b,'mbef,fj->bmje');
         




end

