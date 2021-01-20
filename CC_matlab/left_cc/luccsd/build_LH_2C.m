function [X_abij] = build_LH_2C(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_jacobi,flag_ground)
    
    H1B = HBar_t.H1B;
    H2B = HBar_t.H2B; H2C = HBar_t.H2C;
 %   H3C = HBar_t.H3C; H3D = HBar_t.H3D;
    
    t1a = cc_t.t1a; t1b = cc_t.t1b;
    t2a = cc_t.t2a; t2b = cc_t.t2b; t2c = cc_t.t2c;
    
    % 
    Nocc_a = sys.Nocc_alpha; Nocc_b = sys.Nocc_alpha;
    Nunocc_a = sys.Nvir_alpha; Nunocc_b = sys.Nvir_beta;
    Zunocc_b = ones(Nunocc_b) - eye(Nunocc_b);
    Zocc_b = ones(Nocc_b) - eye(Nocc_b);
    
    % permute from abij into ijab
    l1b = permute(l1b,[2,1]);
    l2b = permute(l2b,[3,4,1,2]); l2c = permute(l2c,[3,4,1,2]);
    
    X_ijab = zeros(Nocc_b,Nocc_b,Nunocc_b,Nunocc_b);

    if flag_jacobi
        X_ijab = X_ijab + einsum_kg(H1B.vv.*Zunocc_b,l2c,'ea,ijeb->ijab')...
                        - einsum_kg(H1B.vv.*Zunocc_b,l2c,'eb,ijea->ijab');

        X_ijab = X_ijab - einsum_kg(H1B.oo.*Zocc_b,l2c,'im,mjab->ijab')...
                        + einsum_kg(H1B.oo.*Zocc_b,l2c,'jm,miab->ijab');
%         X_ijab = X_ijab + einsum_kg(H1B.vv,l2c,'ea,ijeb->ijab')...
%                         - einsum_kg(H1B.vv,l2c,'eb,ijea->ijab');
% 
%         X_ijab = X_ijab - einsum_kg(H1B.oo,l2c,'im,mjab->ijab')...
%                         + einsum_kg(H1B.oo,l2c,'jm,miab->ijab');
    else
        X_ijab = X_ijab + einsum_kg(H1B.vv,l2c,'ea,ijeb->ijab')...
                        - einsum_kg(H1B.vv,l2c,'eb,ijea->ijab');

        X_ijab = X_ijab - einsum_kg(H1B.oo,l2c,'im,mjab->ijab')...
                        + einsum_kg(H1B.oo,l2c,'jm,miab->ijab');
    end
                
    X_ijab = X_ijab - einsum_kg(H2C.ooov,l1b,'ijmb,ma->ijab')...
                    + einsum_kg(H2C.ooov,l1b,'ijma,mb->ijab');
               
    X_ijab = X_ijab + einsum_kg(H2C.vovv,l1b,'ejab,ie->ijab')...
                    - einsum_kg(H2C.vovv,l1b,'eiab,je->ijab');
%                 
    X_ijab = X_ijab + 0.5*einsum_kg(H2C.vvvv,l2c,'efab,ijef->ijab');    
    X_ijab = X_ijab + 0.5*einsum_kg(H2C.oooo,l2c,'ijmn,mnab->ijab');
%     
    X_ijab = X_ijab + einsum_kg(H2C.voov,l2c,'ejmb,imae->ijab')...
                    - einsum_kg(H2C.voov,l2c,'eimb,jmae->ijab')...
                    - einsum_kg(H2C.voov,l2c,'ejma,imbe->ijab')...
                    + einsum_kg(H2C.voov,l2c,'eima,jmbe->ijab');
                
    X_ijab = X_ijab + einsum_kg(H2B.voov,l2b,'ejmb,miea->ijab')...
                    - einsum_kg(H2B.voov,l2b,'eimb,mjea->ijab')...
                    - einsum_kg(H2B.voov,l2b,'ejma,mieb->ijab')...
                    + einsum_kg(H2B.voov,l2b,'eima,mjeb->ijab');
                
    X_ijab = X_ijab - einsum_kg(einsum_kg(l2b,t2b,'nmfa,fenm->ea'),H2C.oovv,'ea,ijeb->ijab')...
                    + einsum_kg(einsum_kg(l2b,t2b,'nmfb,fenm->eb'),H2C.oovv,'eb,ijea->ijab');
               
    X_ijab = X_ijab - 0.5*einsum_kg(einsum_kg(l2c,t2c,'mnaf,efmn->ea'),H2C.oovv,'ea,ijeb->ijab')...
                    + 0.5*einsum_kg(einsum_kg(l2c,t2c,'mnbf,efmn->eb'),H2C.oovv,'eb,ijea->ijab');
                
    X_ijab = X_ijab - einsum_kg(einsum_kg(l2b,t2b,'nife,fenm->im'),H2C.oovv,'im,mjab->ijab')...
                    + einsum_kg(einsum_kg(l2b,t2b,'njfe,fenm->jm'),H2C.oovv,'jm,miab->ijab');
                
    X_ijab = X_ijab - 0.5*einsum_kg(einsum_kg(l2c,t2c,'inef,efmn->im'),H2C.oovv,'im,mjab->ijab')...
                    + 0.5*einsum_kg(einsum_kg(l2c,t2c,'jnef,efmn->jm'),H2C.oovv,'jm,miab->ijab');
%                 
%     X_ijab = X_ijab - 0.5*einsum_kg(H3D.oovovo,l2c,'ijfmbn,mnaf->ijab')...
%                     + 0.5*einsum_kg(H3D.oovovo,l2c,'ijfman,mnbf->ijab');
%                 
%     X_ijab = X_ijab - einsum_kg(H3C.voooov,l2b,'fijnmb,nmfa->ijab')...
%                     + einsum_kg(H3C.voooov,l2b,'fijnma,nmfb->ijab');
%                 
%     X_ijab = X_ijab + 0.5*einsum_kg(H3D.vvovov,l2c,'efjanb,inef->ijab')...
%                     - 0.5*einsum_kg(H3D.vvovov,l2c,'efianb,jnef->ijab');
%                 
%     X_ijab = X_ijab + einsum_kg(H3C.vvoovv,l2b,'fejnab,nife->ijab')...
%                     - einsum_kg(H3C.vvoovv,l2b,'feinab,njfe->ijab');
                
    X_ijab = X_ijab + einsum_kg(H1B.ov,l1b,'jb,ia->ijab')...
                    - einsum_kg(H1B.ov,l1b,'ib,ja->ijab')...
                    - einsum_kg(H1B.ov,l1b,'ja,ib->ijab')...
                    + einsum_kg(H1B.ov,l1b,'ia,jb->ijab');
                    
    % only for ground state
    if flag_ground
        X_ijab = X_ijab + sys.vC_oovv;
    end
    
    X_abij = permute(X_ijab,[3,4,1,2]);
    
%     l2c_new = zeros(Nunocc_b, Nunocc_b, Nocc_b, Nocc_b);
%     for i = 1:Nocc_b
%         for a = 1:Nunocc_b
%             for j = i+1:Nocc_b
%                 for b = a+1:Nunocc_b
%                     
%                     denom = H1B.vv(a,a) - H1B.oo(i,i) + H1B.vv(b,b) - H1B.oo(j,j);
%                       
%                     l2c_new(a,b,i,j) = -X_abij(a,b,i,j)/(denom - omega + shift);                   
%                     l2c_new(b,a,i,j) = -l2c_new(a,b,i,j);
%                     l2c_new(a,b,j,i) = -l2c_new(a,b,i,j);
%                     l2c_new(b,a,j,i) = l2c_new(a,b,i,j);
%                 end
%             end
%         end
%     end
    
    
end