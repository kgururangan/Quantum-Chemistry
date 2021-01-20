function [X_ai] = build_LH_1B(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_jacobi,flag_ground)

    
    H1B = HBar_t.H1B;
    H2B = HBar_t.H2B; H2C = HBar_t.H2C;
%    H3B = HBar_t.H3B; H3C = HBar_t.H3C; H3D = HBar_t.H3D;
        
    t1a = cc_t.t1a; t1b = cc_t.t1b;
    t2a = cc_t.t2a; t2b = cc_t.t2b; t2c = cc_t.t2c;

    % 
    Nocc_b = sys.Nocc_alpha;
    Nunocc_b = sys.Nvir_beta;
    Zunocc_b = ones(Nunocc_b) - eye(Nunocc_b);
    Zocc_b = ones(Nocc_b) - eye(Nocc_b);
    
    % permute from abij into ijab
    l1a = permute(l1a,[2,1]); l1b = permute(l1b,[2,1]);
    l2a = permute(l2a,[3,4,1,2]); l2b = permute(l2b,[3,4,1,2]); l2c = permute(l2c,[3,4,1,2]);
    
    X_ia = zeros(Nocc_b, Nunocc_b);
    if flag_jacobi
        X_ia = X_ia + einsum_kg(H1B.vv.*Zunocc_b,l1b,'ea,ie->ia');
        X_ia = X_ia - einsum_kg(H1B.oo.*Zocc_b,l1b,'im,ma->ia');
%         X_ia = X_ia + einsum_kg(H1B.vv,l1b,'ea,ie->ia');
%         X_ia = X_ia - einsum_kg(H1B.oo,l1b,'im,ma->ia');
    else
        X_ia = X_ia + einsum_kg(H1B.vv,l1b,'ea,ie->ia');
        X_ia = X_ia - einsum_kg(H1B.oo,l1b,'im,ma->ia');
    end
    X_ia = X_ia + einsum_kg(H2B.voov,l1a,'eima,me->ia');
    X_ia = X_ia + einsum_kg(H2C.voov,l1b,'eima,me->ia');
    X_ia = X_ia - 0.5*einsum_kg(H2C.ovoo,l2c,'ifmn,mnaf->ia');
    X_ia = X_ia - einsum_kg(H2B.vooo,l2b,'finm,nmfa->ia');
    X_ia = X_ia + einsum_kg(H2B.vvov,l2b,'fena,nife->ia');
    X_ia = X_ia + 0.5*einsum_kg(H2C.vvvo,l2c,'efan,inef->ia');
    
    
%     X_ia = X_ia + 0.25*einsum_kg(H3D.vvooov,l2c,'efimna,mnef->ia');
    X_ia = X_ia + 0.25*einsum_kg(einsum_kg(l2c,t2c,'mnef,fgnm->ge'),H2C.vovv,'ge,eiga->ia')...
                - 0.25*einsum_kg(einsum_kg(l2c,t2c,'mnef,egnm->gf'),H2C.vovv,'gf,figa->ia')...
                - 0.25*einsum_kg(einsum_kg(l2c,t2c,'mnef,efon->mo'),H2C.ooov,'mo,oima->ia')...
                + 0.25*einsum_kg(einsum_kg(l2c,t2c,'mnef,efom->no'),H2C.ooov,'no,oina->ia');
    
%     X_ia = X_ia + 0.25*einsum_kg(H3B.vvooov,l2a,'efimna,mnef->ia');
    X_ia = X_ia + 0.25*einsum_kg(einsum_kg(l2a,t2a,'mnef,fgnm->ge'),H2B.vovv,'ge,eiga->ia')...
                - 0.25*einsum_kg(einsum_kg(l2a,t2a,'mnef,egnm->gf'),H2B.vovv,'gf,figa->ia')...
                - 0.25*einsum_kg(einsum_kg(l2a,t2a,'mnef,efon->mo'),H2B.ooov,'mo,oima->ia')...
                + 0.25*einsum_kg(einsum_kg(l2a,t2a,'mnef,efom->no'),H2B.ooov,'no,oina->ia');
            
%     X_ia = X_ia + einsum_kg(H3C.vvooov,l2b,'efimna,mnef->ia');
    X_ia = X_ia + einsum_kg(einsum_kg(l2b,t2b,'mnef,gfmn->ge'),H2B.vovv,'ge,eiga->ia')...
                + einsum_kg(einsum_kg(l2b,t2b,'nmfe,fgnm->ge'),H2C.vovv,'ge,eiga->ia')...
                - einsum_kg(einsum_kg(l2b,t2b,'mnef,efon->mo'),H2B.ooov,'mo,oima->ia')...
                - einsum_kg(einsum_kg(l2b,t2b,'nmfe,feno->mo'),H2C.ooov,'mo,oima->ia');

    % only for ground state
    if flag_ground
        X_ia = X_ia + H1B.ov;
    end
    
    X_ai = permute(X_ia,[2,1]);
            
%     l1b_new = zeros(Nunocc_b,Nocc_b);
%     for i = 1:Nocc_b
%         for a = 1:Nunocc_b
%             denom = H1B.vv(a,a) - H1B.oo(i,i);
%             l1b_new(a,i) = -X_ai(a,i)/(denom - omega + shift); 
%         end
%     end

end

