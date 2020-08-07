function [X_abij] = build_LH_2A(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_jacobi,flag_ground)
    
    
    H1A = HBar_t.H1A; 
    H2A = HBar_t.H2A; H2B = HBar_t.H2B; 
%    H3A = HBar_t.H3A; H3B = HBar_t.H3B; 
    
    t1a = cc_t.t1a; t1b = cc_t.t1b;
    t2a = cc_t.t2a; t2b = cc_t.t2b; t2c = cc_t.t2c;
    
    % 
    Nocc_a = sys.Nocc_alpha; 
    Nunocc_a = sys.Nvir_alpha; 
    Zunocc_a = ones(Nunocc_a) - eye(Nunocc_a);
    Zocc_a = ones(Nocc_a) - eye(Nocc_a);
    
    % permute from abij into ijab
    l1a = permute(l1a,[2,1]); 
    l2a = permute(l2a,[3,4,1,2]); l2b = permute(l2b,[3,4,1,2]); 
    
    X_ijab = zeros(Nocc_a,Nocc_a,Nunocc_a,Nunocc_a);
    if flag_jacobi
        X_ijab = X_ijab + einsum_kg(H1A.vv.*Zunocc_a,l2a,'ea,ijeb->ijab')...
                        - einsum_kg(H1A.vv.*Zunocc_a,l2a,'eb,ijea->ijab');
        X_ijab = X_ijab - einsum_kg(H1A.oo.*Zocc_a,l2a,'im,mjab->ijab')...
                        + einsum_kg(H1A.oo.*Zocc_a,l2a,'jm,miab->ijab');
    else
        X_ijab = X_ijab + einsum_kg(H1A.vv,l2a,'ea,ijeb->ijab')...
                        - einsum_kg(H1A.vv,l2a,'eb,ijea->ijab');
        X_ijab = X_ijab - einsum_kg(H1A.oo,l2a,'im,mjab->ijab')...
                        + einsum_kg(H1A.oo,l2a,'jm,miab->ijab');
    end
    X_ijab = X_ijab + einsum_kg(H1A.ov,l1a,'jb,ia->ijab')...
                    - einsum_kg(H1A.ov,l1a,'ja,ib->ijab')...
                    - einsum_kg(H1A.ov,l1a,'ib,ja->ijab')...
                    + einsum_kg(H1A.ov,l1a,'ia,jb->ijab');
                
    X_ijab = X_ijab - 0.5*einsum_kg(einsum_kg(l2a,t2a,'mnaf,efmn->ea'),H2A.oovv,'ea,ijeb->ijab') ...
                    + 0.5*einsum_kg(einsum_kg(l2a,t2a,'mnbf,efmn->eb'),H2A.oovv,'eb,ijea->ijab');
                
    X_ijab = X_ijab - einsum_kg(einsum_kg(l2b,t2b,'mnaf,efmn->ea'),H2A.oovv,'ea,ijeb->ijab')...
                    + einsum_kg(einsum_kg(l2b,t2b,'mnbf,efmn->eb'),H2A.oovv,'eb,ijea->ijab');
                
    X_ijab = X_ijab - 0.5*einsum_kg(einsum_kg(l2a,t2a,'inef,efmn->im'),H2A.oovv,'im,mjab->ijab')...
                    + 0.5*einsum_kg(einsum_kg(l2a,t2a,'jnef,efmn->jm'),H2A.oovv,'jm,miab->ijab');
                
    X_ijab = X_ijab - einsum_kg(einsum_kg(l2b,t2b,'inef,efmn->im'),H2A.oovv,'im,mjab->ijab')...
                    + einsum_kg(einsum_kg(l2b,t2b,'jnef,efmn->jm'),H2A.oovv,'jm,miab->ijab');
                
%     X_ijab = X_ijab - 0.5*einsum_kg(einsum_kg(sys.vA_oovv,t2a,'mnfe,afij->amnije'),l2a,'fijnmb,mnaf->ijab')...
%                     + 0.5*einsum_kg(einsum_kg(sys.vA_oovv,t2a,'mnfe,afij->amnije'),l2a,'fijnma,mnbf->ijab');
%                 
%     X_ijab = X_ijab - einsum_kg(einsum_kg(sys.vA_oovv,t2b,'mnfe,faji->mnajei'),l2b,'ijfmbn,mnaf->ijab')...
%                     + einsum_kg(einsum_kg(sys.vA_oovv,t2b,'mnfe,faji->mnajei'),l2b,'ijfman,mnbf->ijab');
% 
% %     X_ijab = X_ijab + 0.5*einsum_kg(H3A.vvovov,l2a,'efjanb,inef->ijab')...
% %                     - 0.5*einsum_kg(H3A.vvovov,l2a,'efianb,jnef->ijab');
% %                 
% %     X_ijab = X_ijab + einsum_kg(H3B.vovvvo,l2b,'ejfabn,inef->ijab')...
% %                     - einsum_kg(H3B.vovvvo,l2b,'eifabn,jnef->ijab');
%                 
%     X_ijab = X_ijab + 0.5*einsum_kg(H3A.vvovov,l2a,'efjanb,inef->ijab')...
%                     - 0.5*einsum_kg(H3A.vvovov,l2a,'efianb,jnef->ijab');
%                 
%     X_ijab = X_ijab + einsum_kg(H3B.vovvvo,l2b,'ejfabn,inef->ijab')...
%                     - einsum_kg(H3B.vovvvo,l2b,'eifabn,jnef->ijab');

    % THIS WORKS WHEN H2A.ovvo = -permute(H2A.voov,[2,1,4,3]): 
    % WHY - SIGN???
    X_ijab = X_ijab + einsum_kg(H2A.voov,l2a,'eima,mjeb->ijab')...
                    - einsum_kg(H2A.voov,l2a,'ejma,mieb->ijab')...
                    - einsum_kg(H2A.voov,l2a,'eimb,mjea->ijab')...
                    + einsum_kg(H2A.voov,l2a,'ejmb,miea->ijab');
%     X_ijab = X_ijab + einsum_kg(H2A.ovvo,l2a,'ieam,mjeb->ijab')...
%                     - einsum_kg(H2A.ovvo,l2a,'jeam,mieb->ijab')...
%                     - einsum_kg(H2A.ovvo,l2a,'iebm,mjea->ijab')...
%                     + einsum_kg(H2A.ovvo,l2a,'jebm,miea->ijab');
                
    X_ijab = X_ijab + einsum_kg(H2B.ovvo,l2b,'ieam,jmbe->ijab')...
                    - einsum_kg(H2B.ovvo,l2b,'jeam,imbe->ijab')...
                    - einsum_kg(H2B.ovvo,l2b,'iebm,jmae->ijab')...
                    + einsum_kg(H2B.ovvo,l2b,'jebm,imae->ijab');
                
    X_ijab = X_ijab + 0.5*einsum_kg(H2A.oooo,l2a,'ijmn,mnab->ijab');
    X_ijab = X_ijab + 0.5*einsum_kg(H2A.vvvv,l2a,'efab,ijef->ijab');
    X_ijab = X_ijab + einsum_kg(H2A.vovv,l1a,'ejab,ie->ijab')...
                    - einsum_kg(H2A.vovv,l1a,'eiab,je->ijab');
    X_ijab = X_ijab - einsum_kg(H2A.ooov,l1a,'ijmb,ma->ijab')...
                    + einsum_kg(H2A.ooov,l1a,'ijma,mb->ijab');
%                 
    % only for ground state
    if flag_ground
        X_ijab = X_ijab + sys.vA_oovv;
    end
    
    X_abij = permute(X_ijab,[3,4,1,2]);
    
%     l2a_new = zeros(Nunocc_a, Nunocc_a, Nocc_a, Nocc_a);
%     for i = 1:Nocc_a
%         for a = 1:Nunocc_a
%             for j = i+1:Nocc_a
%                 for b = a+1:Nunocc_a
%                     
%                     denom = H1A.vv(a,a) - H1A.oo(i,i) + H1A.vv(b,b) - H1A.oo(j,j);
%                       
%                     l2a_new(a,b,i,j) = -X_abij(a,b,i,j)/(denom  - omega + shift);                   
%                     l2a_new(b,a,i,j) = -l2a_new(a,b,i,j);
%                     l2a_new(a,b,j,i) = -l2a_new(a,b,i,j);
%                     l2a_new(b,a,j,i) = l2a_new(a,b,i,j);
%                 end
%             end
%         end
%     end
%     
    
end