function [L3D] = build_L3D_approx(cc_t,HBar_t,sys,nroot)

    fprintf('Approximate L3D construction...')

    tic
    
    if nroot == 0
        nroot = 1;
    else
        nroot = nroot + 1;
    end
    l1b = cc_t.l1b{nroot}; l2c = cc_t.l2c{nroot};
    
    H2C = HBar_t.H2C; H1B = HBar_t.H1B;
    
    % permute L1 and L2 into ijab
    l1b = permute(l1b,[2,1]); 
    l2c = permute(l2c,[3,4,1,2]);    

%     DC1 =  einsum_kg(l1b,H2C.oovv,'kc,ijab->ijkabc') ... % (1)
%          - einsum_kg(l1b,H2C.oovv,'ka,ijcb->ijkabc') ... % (ac)
%          - einsum_kg(l1b,H2C.oovv,'kb,ijac->ijkabc') ... % (bc)
%          - einsum_kg(l1b,H2C.oovv,'ic,kjab->ijkabc') ... % (ik)
%          - einsum_kg(l1b,H2C.oovv,'jc,ikab->ijkabc') ... % (jk)
%          + einsum_kg(l1b,H2C.oovv,'ia,kjcb->ijkabc') ... % (ik)(ac)
%          + einsum_kg(l1b,H2C.oovv,'ib,kjac->ijkabc') ... % (ik)(bc)
%          + einsum_kg(l1b,H2C.oovv,'ja,ikcb->ijkabc') ... % (jk)(ac)
%          + einsum_kg(l1b,H2C.oovv,'jb,ikac->ijkabc');    % (jk)(bc)
    
    DC1 = einsum_kg(l1b,H2C.oovv,'kc,ijab->ijkabc');
    DC1 = DC1 - permute(DC1,[1,3,2,4,5,6]) - permute(DC1,[3,2,1,4,5,6])...
              - permute(DC1,[1,2,3,6,5,4]) - permute(DC1,[1,2,3,4,6,5])...
              + permute(DC1,[1,3,2,6,5,4]) + permute(DC1,[1,3,2,4,6,5])...
              + permute(DC1,[3,2,1,6,5,4]) + permute(DC1,[3,2,1,4,6,5]);

%     DC2 =  einsum_kg(H1B.ov,l2c,'kc,ijab->ijkabc') ... % (1)
%          - einsum_kg(H1B.ov,l2c,'ka,ijcb->ijkabc') ... % (ac)
%          - einsum_kg(H1B.ov,l2c,'kb,ijac->ijkabc') ... % (bc)
%          - einsum_kg(H1B.ov,l2c,'ic,kjab->ijkabc') ... % (ik)
%          - einsum_kg(H1B.ov,l2c,'jc,ikab->ijkabc') ... % (jk)
%          + einsum_kg(H1B.ov,l2c,'ia,kjcb->ijkabc') ... % (ik)(ac)
%          + einsum_kg(H1B.ov,l2c,'ib,kjac->ijkabc') ... % (ik)(bc)
%          + einsum_kg(H1B.ov,l2c,'ja,ikcb->ijkabc') ... % (jk)(ac)
%          + einsum_kg(H1B.ov,l2c,'jb,ikac->ijkabc');    % (jk)(bc)

    DC2 = einsum_kg(H1B.ov,l2c,'kc,ijab->ijkabc');
    DC2 = DC2 - permute(DC2,[1,3,2,4,5,6]) - permute(DC2,[3,2,1,4,5,6])...
              - permute(DC2,[1,2,3,6,5,4]) - permute(DC2,[1,2,3,4,6,5])...
              + permute(DC2,[1,3,2,6,5,4]) + permute(DC2,[1,3,2,4,6,5])...
              + permute(DC2,[3,2,1,6,5,4]) + permute(DC2,[3,2,1,4,6,5]);


%     D1 =   einsum_kg(H2C.ovvv,l2c,'ieab,jkec->ijkabc') ... % (1)
%          - einsum_kg(H2C.ovvv,l2c,'jeab,ikec->ijkabc') ... % (ij)
%          - einsum_kg(H2C.ovvv,l2c,'keab,jiec->ijkabc') ... % (ik)
%          - einsum_kg(H2C.ovvv,l2c,'iecb,jkea->ijkabc') ... % (ac)
%          - einsum_kg(H2C.ovvv,l2c,'ieac,jkeb->ijkabc') ... % (bc)
%          + einsum_kg(H2C.ovvv,l2c,'jecb,ikea->ijkabc') ... % (ij)(ac)
%          + einsum_kg(H2C.ovvv,l2c,'jeac,ikeb->ijkabc') ... % (ij)(bc)
%          + einsum_kg(H2C.ovvv,l2c,'kecb,jiea->ijkabc') ... % (ik)(ac)
%          + einsum_kg(H2C.ovvv,l2c,'keac,jieb->ijkabc');    % (ik)(bc)
    D1 = einsum_kg(H2C.ovvv,l2c,'ieab,jkec->ijkabc');
    D1 = D1   - permute(D1,[2,1,3,4,5,6]) - permute(D1,[3,2,1,4,5,6])...
              - permute(D1,[1,2,3,6,5,4]) - permute(D1,[1,2,3,4,6,5])...
              + permute(D1,[2,1,3,6,5,4]) + permute(D1,[2,1,3,4,6,5])...
              + permute(D1,[3,2,1,6,5,4]) + permute(D1,[3,2,1,4,6,5]);

%     D2 =   einsum_kg(H2C.oovo,l2c,'ijam,mkbc->ijkabc') ... % (1)
%          - einsum_kg(H2C.oovo,l2c,'kjam,mibc->ijkabc') ... % (ik)
%          - einsum_kg(H2C.oovo,l2c,'ikam,mjbc->ijkabc') ... % (jk)
%          - einsum_kg(H2C.oovo,l2c,'ijbm,mkac->ijkabc') ... % (ab)
%          - einsum_kg(H2C.oovo,l2c,'ijcm,mkba->ijkabc') ... % (ac)
%          + einsum_kg(H2C.oovo,l2c,'kjbm,miac->ijkabc') ... % (ik)(ab)
%          + einsum_kg(H2C.oovo,l2c,'kjcm,miba->ijkabc') ... % (ik)(ac)
%          + einsum_kg(H2C.oovo,l2c,'ikbm,mjac->ijkabc') ... % (jk)(ab)
%          + einsum_kg(H2C.oovo,l2c,'ikcm,mjba->ijkabc');    % (jk)(ac)
    D2 = einsum_kg(H2C.oovo,l2c,'ijam,mkbc->ijkabc');
    D2 = D2 - permute(D2,[3,2,1,4,5,6]) - permute(D2,[1,3,2,4,5,6]) ...
            - permute(D2,[1,2,3,5,4,6]) - permute(D2,[1,2,3,6,5,4]) ...
            + permute(D2,[3,2,1,5,4,6]) + permute(D2,[3,2,1,6,5,4]) ...
            + permute(D2,[1,3,2,5,4,6]) + permute(D2,[1,3,2,6,5,4]);

    L3D = DC1 + DC2 + D1 - D2;

    fprintf(' finished in %4.2f s\n',toc)
end

