function [X_ai] = build_LH_1A(l1a,l1b,l2a,l2b,l2c,HBar_t,cc_t,sys,flag_jacobi,flag_ground)

    
    H1A = HBar_t.H1A; 
    H2A = HBar_t.H2A; H2B = HBar_t.H2B; 
%    H3A = HBar_t.H3A; H3B = HBar_t.H3B; H3C = HBar_t.H3C; 
    
    t1a = cc_t.t1a; t1b = cc_t.t1b;
    t2a = cc_t.t2a; t2b = cc_t.t2b; t2c = cc_t.t2c;
    
    % 
    Nocc_a = sys.Nocc_alpha; 
    Nunocc_a = sys.Nvir_alpha; 
    Zunocc_a = ones(Nunocc_a) - eye(Nunocc_a);
    Zocc_a = ones(Nocc_a) - eye(Nocc_a);
    
    % permute from ai into ia and abij into ijab
    l1a = permute(l1a,[2,1]); l1b = permute(l1b,[2,1]);
    l2a = permute(l2a,[3,4,1,2]); l2b = permute(l2b,[3,4,1,2]); l2c = permute(l2c,[3,4,1,2]);
    
    X_ia = zeros(Nocc_a, Nunocc_a);
    if flag_jacobi
        X_ia = X_ia + einsum_kg(H1A.vv.*Zunocc_a,l1a,'ea,ie->ia');
        X_ia = X_ia - einsum_kg(H1A.oo.*Zocc_a,l1a,'im,ma->ia');
    else
        X_ia = X_ia + einsum_kg(H1A.vv,l1a,'ea,ie->ia');
        X_ia = X_ia - einsum_kg(H1A.oo,l1a,'im,ma->ia');
    end
    X_ia = X_ia + einsum_kg(H2A.voov,l1a,'eima,me->ia');
    X_ia = X_ia + einsum_kg(H2B.ovvo,l1b,'ieam,me->ia');
    X_ia = X_ia + 0.5*einsum_kg(H2A.vvvo,l2a,'efan,inef->ia');
    X_ia = X_ia + einsum_kg(H2B.vvvo,l2b,'efan,inef->ia');
    X_ia = X_ia - 0.5*einsum_kg(H2A.ovoo,l2a,'ifmn,mnaf->ia');
    X_ia = X_ia - einsum_kg(H2B.ovoo,l2b,'ifmn,mnaf->ia');


%     X_ia = X_ia + 0.25*einsum_kg(H3A.ovvvoo,l2a,'iefamn,mnef->ia');
    
    X_ia = X_ia + 0.25*einsum_kg(einsum_kg(l2a,t2a,'mnef,fgnm->ge'),H2A.vovv,'ge,eiga->ia') ...
                - 0.25*einsum_kg(einsum_kg(l2a,t2a,'mnef,egnm->gf'),H2A.vovv,'gf,figa->ia') ...
                - 0.25*einsum_kg(einsum_kg(l2a,t2a,'moef,efno->mn'),H2A.ooov,'mn,nima->ia') ...
                + 0.25*einsum_kg(einsum_kg(l2a,t2a,'moef,efnm->on'),H2A.ooov,'on,nioa->ia');
%     
%      X_ia = X_ia + einsum_kg(H3B.ovvvoo,l2b,'iefamn,mnef->ia');

    X_ia = X_ia - einsum_kg(einsum_kg(l2b,t2b,'ijab,abin->jn'),H2B.oovo,'jn,mnej->me') ...
                + einsum_kg(einsum_kg(l2b,t2b,'ijab,afij->fb'),H2B.ovvv,'fb,mbef->me') ...
                + einsum_kg(einsum_kg(l2b,t2b,'ijab,fbij->fa'),H2A.ovvv,'fa,maef->me') ...
                - einsum_kg(einsum_kg(l2b,t2b,'ijab,abnj->in'),H2A.oovo,'in,mnei->me');
%     
%     X_ia = X_ia + 0.25*einsum_kg(H3C.ovvvoo,l2c,'iefamn,mnef->ia');

    X_ia = X_ia + 0.25*einsum_kg(einsum_kg(l2c,t2c,'ijab,fbij->fa'),H2B.ovvv,'fa,maef->me') ...
                - 0.25*einsum_kg(einsum_kg(l2c,t2c,'ijab,faij->fb'),H2B.ovvv,'fb,mbef->me') ...
                - 0.25*einsum_kg(einsum_kg(l2c,t2c,'ijab,abnj->in'),H2B.oovo,'in,mnei->me') ...
                + 0.25*einsum_kg(einsum_kg(l2c,t2c,'ijab,abni->jn'),H2B.oovo,'jn,mnej->me');
    
    % only for ground state
    if flag_ground
        X_ia = X_ia + H1A.ov;
    end

    %X_ia = D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + D9 + D10 + D11 + D12;
    X_ai = permute(X_ia,[2,1]);
    
%     l1a_new = zeros(Nunocc_a, Nocc_a);
%     for i = 1:Nocc_a
%         for a = 1:Nunocc_a
%             denom = H1A.vv(a,a) - H1A.oo(i,i);
%             l1a_new(a,i) = -X_ai(a,i)/(denom - omega + shift); 
%         end
%     end
    
    
end

