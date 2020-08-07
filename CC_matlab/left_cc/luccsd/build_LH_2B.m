function [X_abij] = build_LH_2B(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_jacobi,flag_ground)

    
    H1A = HBar_t.H1A; H1B = HBar_t.H1B;
    H2A = HBar_t.H2A; H2B = HBar_t.H2B; H2C = HBar_t.H2C;
 %   H3B = HBar_t.H3B; H3C = HBar_t.H3C; 
    
    t1a = cc_t.t1a; t1b = cc_t.t1b;
    t2a = cc_t.t2a; t2b = cc_t.t2b; t2c = cc_t.t2c;
    
    % 
    Nocc_a = sys.Nocc_alpha; Nocc_b = sys.Nocc_alpha;
    Nunocc_a = sys.Nvir_alpha; Nunocc_b = sys.Nvir_beta;
    Zunocc_a = ones(Nunocc_a) - eye(Nunocc_a);
    Zunocc_b = ones(Nunocc_b) - eye(Nunocc_b);
    Zocc_a = ones(Nocc_a) - eye(Nocc_a);
    Zocc_b = ones(Nocc_b) - eye(Nocc_b);
    
    % permute from abij into ijab
    l1a = permute(l1a,[2,1]); l1b = permute(l1b,[2,1]);
    l2a = permute(l2a,[3,4,1,2]); l2b = permute(l2b,[3,4,1,2]); l2c = permute(l2c,[3,4,1,2]);
    
    X_ijab = zeros(Nocc_a,Nocc_b,Nunocc_a,Nunocc_b);
    
    X_ijab = X_ijab - einsum_kg(H2B.ooov,l1a,'ijmb,ma->ijab');
    X_ijab = X_ijab - einsum_kg(H2B.oovo,l1b,'ijam,mb->ijab');
    
    X_ijab = X_ijab + einsum_kg(H2B.vovv,l1a,'ejab,ie->ijab');
    X_ijab = X_ijab + einsum_kg(H2B.ovvv,l1b,'ieab,je->ijab');
%     
    X_ijab = X_ijab + einsum_kg(H2B.oooo,l2b,'ijmn,mnab->ijab');
    X_ijab = X_ijab + einsum_kg(H2B.vvvv,l2b,'efab,ijef->ijab');
%     
    X_ijab = X_ijab + einsum_kg(H2B.voov,l2a,'ejmb,imae->ijab');
    X_ijab = X_ijab + einsum_kg(H2A.voov,l2b,'eima,mjeb->ijab');
    X_ijab = X_ijab + einsum_kg(H2C.voov,l2b,'ejmb,imae->ijab');
    X_ijab = X_ijab + einsum_kg(H2B.ovvo,l2c,'ieam,mjeb->ijab');
    X_ijab = X_ijab - einsum_kg(H2B.ovov,l2b,'iemb,mjae->ijab');
    X_ijab = X_ijab - einsum_kg(H2B.vovo,l2b,'ejam,imeb->ijab');
%     
    %X_ijab = X_ijab - 0.5*einsum_kg(H3B.voooov,l2a,'fijnmb,mnaf->ijab');
    X_ijab = X_ijab - 0.5*einsum_kg(einsum_kg(l2a,t2a,'ijab,fbij->fa'),H2B.oovv,'fa,nmfe->nmae');
    %X_ijab = X_ijab - einsum_kg(H3C.ovooov,l2b,'ifjmnb,mnaf->ijab');
    X_ijab = X_ijab - einsum_kg(einsum_kg(l2b,t2b,'mnaf,efmn->ea'),H2B.oovv,'ea,ijeb->ijab');
    %X_ijab = X_ijab - einsum_kg(H3B.vooovo,l2b,'fijnam,nmfb->ijab');
    X_ijab = X_ijab - einsum_kg(einsum_kg(l2b,t2b,'nmfb,fenm->eb'),H2B.oovv,'eb,ijae->ijab');
    %X_ijab = X_ijab - 0.5*einsum_kg(H3C.ovovoo,l2c,'ifjanm,mnbf->ijab');
    X_ijab = X_ijab - 0.5*einsum_kg(einsum_kg(l2c,t2c,'mnbf,efmn->eb'),H2B.oovv,'eb,ijae->ijab');
    
    %X_ijab = X_ijab + 0.5*einsum_kg(H3B.vvovov,l2a,'efjanb,inef->ijab');
    X_ijab = X_ijab - 0.5*einsum_kg(einsum_kg(l2a,t2a,'inef,efmn->im'),H2B.oovv,'im,mjab->ijab');
    %X_ijab = X_ijab + einsum_kg(H3C.vvovov,l2b,'efjanb,inef->ijab');
    X_ijab = X_ijab - einsum_kg(einsum_kg(l2b,t2b,'inef,efmn->im'),H2B.oovv,'im,mjab->ijab');
    %X_ijab = X_ijab + einsum_kg(H3B.ovvvov,l2b,'ifeanb,njfe->ijab');
    X_ijab = X_ijab - einsum_kg(einsum_kg(l2b,t2b,'njfe,fenm->jm'),H2B.oovv,'jm,imab->ijab');
    %X_ijab = X_ijab + 0.5*einsum_kg(H3C.ovvvov,l2c,'ifeanb,njfe->ijab');
    X_ijab = X_ijab - 0.5*einsum_kg(einsum_kg(l2c,t2c,'jnef,efmn->jm'),H2B.oovv,'jm,imab->ijab');
    
    if flag_jacobi
        X_ijab = X_ijab + einsum_kg(H1A.vv.*Zunocc_a,l2b,'ea,ijeb->ijab');
        X_ijab = X_ijab + einsum_kg(H1B.vv.*Zunocc_b,l2b,'eb,ijae->ijab');

        X_ijab = X_ijab - einsum_kg(H1A.oo.*Zocc_a,l2b,'im,mjab->ijab');
        X_ijab = X_ijab - einsum_kg(H1B.oo.*Zocc_b,l2b,'jm,imab->ijab');
    else
        X_ijab = X_ijab + einsum_kg(H1A.vv,l2b,'ea,ijeb->ijab');
        X_ijab = X_ijab + einsum_kg(H1B.vv,l2b,'eb,ijae->ijab');

        X_ijab = X_ijab - einsum_kg(H1A.oo,l2b,'im,mjab->ijab');
        X_ijab = X_ijab - einsum_kg(H1B.oo,l2b,'jm,imab->ijab');
    end
    X_ijab = X_ijab + einsum_kg(H1B.ov,l1a,'jb,ia->ijab');
    X_ijab = X_ijab + einsum_kg(H1A.ov,l1b,'ia,jb->ijab');
              
    % only for ground state
    if flag_ground
        X_ijab = X_ijab + sys.vB_oovv;
    end
    
    X_abij = permute(X_ijab,[3,4,1,2]);
    
%     l2b_new = zeros(Nunocc_a, Nunocc_b, Nocc_a, Nocc_b);
%     for i = 1:Nocc_a
%         for a = 1:Nunocc_a
%             for j = 1:Nocc_b
%                 for b = 1:Nunocc_b                    
%                     denom = H1A.vv(a,a) - H1A.oo(i,i) + H1B.vv(b,b) - H1B.oo(j,j);                      
%                     l2b_new(a,b,i,j) = -X_abij(a,b,i,j)/(denom - omega + shift);                   
%                 end
%             end
%         end
%     end
%     
    
end