function [L3A] = build_L3A_approx(cc_t,HBar_t,sys,nroot)

    fprintf('Approximate L3A construction...')

    tic
    
    if nroot == 0
        l1a = cc_t.l1a{1}; l2a = cc_t.l2a{1};
    else
        l1a = cc_t.l1a{nroot+1}; l2a = cc_t.l2a{nroot+1};
    end
    
    H2A = HBar_t.H2A; H1A = HBar_t.H1A;
    
    % permute L1 and L2 into ijab
    l1a = permute(l1a,[2,1]); 
    l2a = permute(l2a,[3,4,1,2]); 
    

    DC1 =  einsum_kg(l1a,H2A.oovv,'kc,ijab->ijkabc') ... % (1)
         - einsum_kg(l1a,H2A.oovv,'ka,ijcb->ijkabc') ... % (ac)
         - einsum_kg(l1a,H2A.oovv,'kb,ijac->ijkabc') ... % (bc)
         - einsum_kg(l1a,H2A.oovv,'ic,kjab->ijkabc') ... % (ik)
         - einsum_kg(l1a,H2A.oovv,'jc,ikab->ijkabc') ... % (jk)
         + einsum_kg(l1a,H2A.oovv,'ia,kjcb->ijkabc') ... % (ik)(ac)
         + einsum_kg(l1a,H2A.oovv,'ib,kjac->ijkabc') ... % (ik)(bc)
         + einsum_kg(l1a,H2A.oovv,'ja,ikcb->ijkabc') ... % (jk)(ac)
         + einsum_kg(l1a,H2A.oovv,'jb,ikac->ijkabc');    % (jk)(bc)

    DC2 =  einsum_kg(H1A.ov,l2a,'kc,ijab->ijkabc') ... % (1)
         - einsum_kg(H1A.ov,l2a,'ka,ijcb->ijkabc') ... % (ac)
         - einsum_kg(H1A.ov,l2a,'kb,ijac->ijkabc') ... % (bc)
         - einsum_kg(H1A.ov,l2a,'ic,kjab->ijkabc') ... % (ik)
         - einsum_kg(H1A.ov,l2a,'jc,ikab->ijkabc') ... % (jk)
         + einsum_kg(H1A.ov,l2a,'ia,kjcb->ijkabc') ... % (ik)(ac)
         + einsum_kg(H1A.ov,l2a,'ib,kjac->ijkabc') ... % (ik)(bc)
         + einsum_kg(H1A.ov,l2a,'ja,ikcb->ijkabc') ... % (jk)(ac)
         + einsum_kg(H1A.ov,l2a,'jb,ikac->ijkabc');    % (jk)(bc)


    D1 =   einsum_kg(H2A.ovvv,l2a,'ieab,jkec->ijkabc') ... % (1)
         - einsum_kg(H2A.ovvv,l2a,'jeab,ikec->ijkabc') ... % (ij)
         - einsum_kg(H2A.ovvv,l2a,'keab,jiec->ijkabc') ... % (ik)
         - einsum_kg(H2A.ovvv,l2a,'iecb,jkea->ijkabc') ... % (ac)
         - einsum_kg(H2A.ovvv,l2a,'ieac,jkeb->ijkabc') ... % (bc)
         + einsum_kg(H2A.ovvv,l2a,'jecb,ikea->ijkabc') ... % (ij)(ac)
         + einsum_kg(H2A.ovvv,l2a,'jeac,ikeb->ijkabc') ... % (ij)(bc)
         + einsum_kg(H2A.ovvv,l2a,'kecb,jiea->ijkabc') ... % (ik)(ac)
         + einsum_kg(H2A.ovvv,l2a,'keac,jieb->ijkabc');    % (ik)(bc)

    D2 =   einsum_kg(H2A.oovo,l2a,'ijam,mkbc->ijkabc') ... % (1)
         - einsum_kg(H2A.oovo,l2a,'kjam,mibc->ijkabc') ... % (ik)
         - einsum_kg(H2A.oovo,l2a,'ikam,mjbc->ijkabc') ... % (jk)
         - einsum_kg(H2A.oovo,l2a,'ijbm,mkac->ijkabc') ... % (ab)
         - einsum_kg(H2A.oovo,l2a,'ijcm,mkba->ijkabc') ... % (ac)
         + einsum_kg(H2A.oovo,l2a,'kjbm,miac->ijkabc') ... % (ik)(ab)
         + einsum_kg(H2A.oovo,l2a,'kjcm,miba->ijkabc') ... % (ik)(ac)
         + einsum_kg(H2A.oovo,l2a,'ikbm,mjac->ijkabc') ... % (jk)(ab)
         + einsum_kg(H2A.oovo,l2a,'ikcm,mjba->ijkabc');    % (jk)(ac)

    L3A = DC1 + DC2 + D1 - D2;

    fprintf(' finished in %4.2f s\n',toc)
end

