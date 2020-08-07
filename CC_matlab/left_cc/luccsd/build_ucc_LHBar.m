function [LH] = build_ucc_LHBar(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_ground)

        flag_jacobi = false;
        X1A = build_LH_1A(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_jacobi,flag_ground);
        X1B = build_LH_1B(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_jacobi,flag_ground);
        X2A = build_LH_2A(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_jacobi,flag_ground);
        X2B = build_LH_2B(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_jacobi,flag_ground);
        X2C = build_LH_2C(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_jacobi,flag_ground);
        LH = cat(1,X1A(:),X1B(:),X2A(:),X2B(:),X2C(:));

%     H1A = HBar_t.H1A; H1B = HBar_t.H1B;
%     H2A = HBar_t.H2A; H2B = HBar_t.H2B; H2C = HBar_t.H2C;
%     H3A = HBar_t.H3A; H3B = HBar_t.H3B; H3C = HBar_t.H3C; H3D = HBar_t.H3D;
%     
%     [Nunocc_a,Nunocc_b,Nocc_a,Nocc_b] = size(l2b);
% 
%     % permute from ai into ia and abij into ijab
%     l1a = permute(l1a,[2,1]); l1b = permute(l1b,[2,1]);
%     l2a = permute(l2a,[3,4,1,2]); l2b = permute(l2b,[3,4,1,2]); l2c = permute(l2c,[3,4,1,2]);
%     
%     % L1A part
%     X1A = zeros(Nocc_a, Nunocc_a);
%     X1A = X1A + einsum_kg(H1A.vv,l1a,'ea,ie->ia');
%     X1A = X1A - einsum_kg(H1A.oo,l1a,'im,ma->ia');
%     X1A = X1A + einsum_kg(H2A.voov,l1a,'eima,me->ia');
%     X1A = X1A + einsum_kg(H2B.ovvo,l1b,'ieam,me->ia');
%     X1A = X1A + 0.5*einsum_kg(H2A.vvvo,l2a,'efan,inef->ia');
%     X1A = X1A + einsum_kg(H2B.vvvo,l2b,'efan,inef->ia');
%     X1A = X1A - 0.5*einsum_kg(H2A.ovoo,l2a,'ifmn,mnaf->ia');
%     X1A = X1A - einsum_kg(H2B.ovoo,l2b,'ifmn,mnaf->ia');
%     X1A = X1A + 0.25*einsum_kg(H3A.ovvvoo,l2a,'iefamn,mnef->ia');
%     X1A = X1A + einsum_kg(H3B.ovvvoo,l2b,'iefamn,mnef->ia');
%     X1A = X1A + 0.25*einsum_kg(H3C.ovvvoo,l2c,'iefamn,mnef->ia');
%     if omega == 0.0
%         X1A = X1A + H1A.ov;
%     end
%     
%     % L1B part
%     X1B = zeros(Nocc_b, Nunocc_b);
%     X1B = X1B + einsum_kg(H1B.vv,l1b,'ea,ie->ia');
%     X1B = X1B - einsum_kg(H1B.oo,l1b,'im,ma->ia');
%     X1B = X1B + einsum_kg(H2B.voov,l1a,'eima,me->ia');
%     X1B = X1B + einsum_kg(H2C.voov,l1b,'eima,me->ia');
%     X1B = X1B - 0.5*einsum_kg(H2C.ovoo,l2c,'ifmn,mnaf->ia');
%     X1B = X1B - einsum_kg(H2B.vooo,l2b,'finm,nmfa->ia');
%     X1B = X1B + einsum_kg(H2B.vvov,l2b,'fena,nife->ia');
%     X1B = X1B + 0.5*einsum_kg(H2C.vvvo,l2c,'efan,inef->ia');
%     X1B = X1B + 0.25*einsum_kg(H3D.vvooov,l2c,'efimna,mnef->ia');
%     X1B = X1B + 0.25*einsum_kg(H3B.vvooov,l2a,'efimna,mnef->ia');
%     X1B = X1B + einsum_kg(H3C.vvooov,l2b,'efimna,mnef->ia');
%     if omega == 0.0
%         X1B = X1B + H1B.ov;
%     end
% 
%     % L2A part
%     X2A = zeros(Nocc_a,Nocc_a,Nunocc_a,Nunocc_a);
%     X2A = X2A + einsum_kg(H1A.vv,l2a,'ea,ijeb->ijab')...
%                     - einsum_kg(H1A.vv,l2a,'eb,ijea->ijab');
%     X2A = X2A - einsum_kg(H1A.oo,l2a,'im,mjab->ijab')...
%                     + einsum_kg(H1A.oo,l2a,'jm,miab->ijab');
%     X2A = X2A + einsum_kg(H1A.ov,l1a,'jb,ia->ijab')...
%                     - einsum_kg(H1A.ov,l1a,'ja,ib->ijab')...
%                     - einsum_kg(H1A.ov,l1a,'ib,ja->ijab')...
%                     + einsum_kg(H1A.ov,l1a,'ia,jb->ijab');               
%     X2A = X2A - 0.5*einsum_kg(H3A.voooov,l2a,'fijnmb,mnaf->ijab')...
%                     + 0.5*einsum_kg(H3A.voooov,l2a,'fijnma,mnbf->ijab');                
%     X2A = X2A - einsum_kg(H3B.oovovo,l2b,'ijfmbn,mnaf->ijab')...
%                     + einsum_kg(H3B.oovovo,l2b,'ijfman,mnbf->ijab');               
%     X2A = X2A + 0.5*einsum_kg(H3A.vvovov,l2a,'efjanb,inef->ijab')...
%                     - 0.5*einsum_kg(H3A.vvovov,l2a,'efianb,jnef->ijab');                
%     X2A = X2A + einsum_kg(H3B.vovvvo,l2b,'ejfabn,inef->ijab')...
%                     - einsum_kg(H3B.vovvvo,l2b,'eifabn,jnef->ijab');
%     X2A = X2A + einsum_kg(H2A.voov,l2a,'eima,mjeb->ijab')...
%                     - einsum_kg(H2A.voov,l2a,'ejma,mieb->ijab')...
%                     - einsum_kg(H2A.voov,l2a,'eimb,mjea->ijab')...
%                     + einsum_kg(H2A.voov,l2a,'ejmb,miea->ijab');             
%     X2A = X2A + einsum_kg(H2B.ovvo,l2b,'ieam,jmbe->ijab')...
%                     - einsum_kg(H2B.ovvo,l2b,'jeam,imbe->ijab')...
%                     - einsum_kg(H2B.ovvo,l2b,'iebm,jmae->ijab')...
%                     + einsum_kg(H2B.ovvo,l2b,'jebm,imae->ijab');               
%     X2A = X2A + 0.5*einsum_kg(H2A.oooo,l2a,'ijmn,mnab->ijab');
%     X2A = X2A + 0.5*einsum_kg(H2A.vvvv,l2a,'efab,ijef->ijab');
%     X2A = X2A + einsum_kg(H2A.vovv,l1a,'ejab,ie->ijab')...
%                     - einsum_kg(H2A.vovv,l1a,'eiab,je->ijab');
%     X2A = X2A - einsum_kg(H2A.ooov,l1a,'ijmb,ma->ijab')...
%                     + einsum_kg(H2A.ooov,l1a,'ijma,mb->ijab');
%     if omega == 0.0
%         X2A = X2A + sys.vA_oovv;
%     end
% 
%     % L2B part
%     X2B = zeros(Nocc_a,Nocc_b,Nunocc_a,Nunocc_b);    
%     X2B = X2B - einsum_kg(H2B.ooov,l1a,'ijmb,ma->ijab');
%     X2B = X2B - einsum_kg(H2B.oovo,l1b,'ijam,mb->ijab');   
%     X2B = X2B + einsum_kg(H2B.vovv,l1a,'ejab,ie->ijab');
%     X2B = X2B + einsum_kg(H2B.ovvv,l1b,'ieab,je->ijab');     
%     X2B = X2B + einsum_kg(H2B.oooo,l2b,'ijmn,mnab->ijab');
%     X2B = X2B + einsum_kg(H2B.vvvv,l2b,'efab,ijef->ijab');     
%     X2B = X2B + einsum_kg(H2B.voov,l2a,'ejmb,imae->ijab');
%     X2B = X2B + einsum_kg(H2A.voov,l2b,'eima,mjeb->ijab');
%     X2B = X2B + einsum_kg(H2C.voov,l2b,'ejmb,imae->ijab');
%     X2B = X2B + einsum_kg(H2B.ovvo,l2c,'ieam,mjeb->ijab');
%     X2B = X2B - einsum_kg(H2B.ovov,l2b,'iemb,mjae->ijab');
%     X2B = X2B - einsum_kg(H2B.vovo,l2b,'ejam,imeb->ijab');     
%     X2B = X2B - 0.5*einsum_kg(H3B.voooov,l2a,'fijnmb,mnaf->ijab');
%     X2B = X2B - einsum_kg(H3C.ovooov,l2b,'ifjmnb,mnaf->ijab');
%     X2B = X2B - einsum_kg(H3B.vooovo,l2b,'fijnam,nmfb->ijab');
%     X2B = X2B - 0.5*einsum_kg(H3C.ovovoo,l2c,'ifjanm,mnbf->ijab');    
%     X2B = X2B + 0.5*einsum_kg(H3B.vvovov,l2a,'efjanb,inef->ijab');
%     X2B = X2B + einsum_kg(H3C.vvovov,l2b,'efjanb,inef->ijab');
%     X2B = X2B + einsum_kg(H3B.ovvvov,l2b,'ifeanb,njfe->ijab');
%     X2B = X2B + 0.5*einsum_kg(H3C.ovvvov,l2c,'ifeanb,njfe->ijab');
%     X2B = X2B + einsum_kg(H1A.vv,l2b,'ea,ijeb->ijab');
%     X2B = X2B + einsum_kg(H1B.vv,l2b,'eb,ijae->ijab');
%     X2B = X2B - einsum_kg(H1A.oo,l2b,'im,mjab->ijab');
%     X2B = X2B - einsum_kg(H1B.oo,l2b,'jm,imab->ijab');
%     X2B = X2B + einsum_kg(H1B.ov,l1a,'jb,ia->ijab');
%     X2B = X2B + einsum_kg(H1A.ov,l1b,'ia,jb->ijab');
%     if omega == 0.0
%         X2B = X2B + sys.vB_oovv;
%     end
% 
%     % L2C part
%     X2C = zeros(Nocc_b,Nocc_b,Nunocc_b,Nunocc_b);
%     X2C = X2C + einsum_kg(H1B.vv,l2c,'ea,ijeb->ijab')...
%                     - einsum_kg(H1B.vv,l2c,'eb,ijea->ijab');
%     X2C = X2C - einsum_kg(H1B.oo,l2c,'im,mjab->ijab')...
%                     + einsum_kg(H1B.oo,l2c,'jm,miab->ijab');
%     X2C = X2C - einsum_kg(H2C.ooov,l1b,'ijmb,ma->ijab')...
%                     + einsum_kg(H2C.ooov,l1b,'ijma,mb->ijab');               
%     X2C = X2C + einsum_kg(H2C.vovv,l1b,'ejab,ie->ijab')...
%                     - einsum_kg(H2C.vovv,l1b,'eiab,je->ijab');                 
%     X2C = X2C + 0.5*einsum_kg(H2C.vvvv,l2c,'efab,ijef->ijab');    
%     X2C = X2C + 0.5*einsum_kg(H2C.oooo,l2c,'ijmn,mnab->ijab');     
%     X2C = X2C + einsum_kg(H2C.voov,l2c,'ejmb,imae->ijab')...
%                     - einsum_kg(H2C.voov,l2c,'eimb,jmae->ijab')...
%                     - einsum_kg(H2C.voov,l2c,'ejma,imbe->ijab')...
%                     + einsum_kg(H2C.voov,l2c,'eima,jmbe->ijab');                
%     X2C = X2C + einsum_kg(H2B.voov,l2b,'ejmb,miea->ijab')...
%                     - einsum_kg(H2B.voov,l2b,'eimb,mjea->ijab')...
%                     - einsum_kg(H2B.voov,l2b,'ejma,mieb->ijab')...
%                     + einsum_kg(H2B.voov,l2b,'eima,mjeb->ijab');                 
%     X2C = X2C - 0.5*einsum_kg(H3D.oovovo,l2c,'ijfmbn,mnaf->ijab')...
%                     + 0.5*einsum_kg(H3D.oovovo,l2c,'ijfman,mnbf->ijab');               
%     X2C = X2C - einsum_kg(H3C.voooov,l2b,'fijnmb,nmfa->ijab')...
%                     + einsum_kg(H3C.voooov,l2b,'fijnma,nmfb->ijab');                
%     X2C = X2C + 0.5*einsum_kg(H3D.vvovov,l2c,'efjanb,inef->ijab')...
%                     - 0.5*einsum_kg(H3D.vvovov,l2c,'efianb,jnef->ijab');                
%     X2C = X2C + einsum_kg(H3C.vvoovv,l2b,'fejnab,nife->ijab')...
%                     - einsum_kg(H3C.vvoovv,l2b,'feinab,njfe->ijab');                
%     X2C = X2C + einsum_kg(H1B.ov,l1b,'jb,ia->ijab')...
%                     - einsum_kg(H1B.ov,l1b,'ib,ja->ijab')...
%                     - einsum_kg(H1B.ov,l1b,'ja,ib->ijab')...
%                     + einsum_kg(H1B.ov,l1b,'ia,jb->ijab');
%     if omega == 0.0
%         X2C = X2C + sys.vC_oovv;
%     end
%      
%     % permute into ai and abij
%     X1A = permute(X1A,[2,1]);
%     X1B = permute(X1B,[2,1]);
%     X2A = permute(X2A,[3,4,1,2]);
%     X2B = permute(X2B,[3,4,1,2]);
%     X2C = permute(X2C,[3,4,1,2]);
%     
%     LH = cat(1,X1A(:),X1B(:),X2A(:),X2B(:),X2C(:));
%     
% %    LAMBDA_resid = LH - omega*LAMBDA;

end

