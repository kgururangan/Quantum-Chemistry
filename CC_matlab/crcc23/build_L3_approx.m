function [L3] = build_L3_approx(L1,L2,HBar)

    fprintf('Approximate L3 construction...')

    tic

    % permute L1 and L2 from abij into ijab
    L1 = permute(L1,[2,1]);     L2 = permute(L2,[3,4,1,2]);

    DC1 =  einsum_kg(L1,HBar{2}{1,1,2,2},'kc,ijab->ijkabc') ... % (1)
         - einsum_kg(L1,HBar{2}{1,1,2,2},'ka,ijcb->ijkabc') ... % (ac)
         - einsum_kg(L1,HBar{2}{1,1,2,2},'kb,ijac->ijkabc') ... % (bc)
         - einsum_kg(L1,HBar{2}{1,1,2,2},'ic,kjab->ijkabc') ... % (ik)
         - einsum_kg(L1,HBar{2}{1,1,2,2},'jc,ikab->ijkabc') ... % (jk)
         + einsum_kg(L1,HBar{2}{1,1,2,2},'ia,kjcb->ijkabc') ... % (ik)(ac)
         + einsum_kg(L1,HBar{2}{1,1,2,2},'ib,kjac->ijkabc') ... % (ik)(bc)
         + einsum_kg(L1,HBar{2}{1,1,2,2},'ja,ikcb->ijkabc') ... % (jk)(ac)
         + einsum_kg(L1,HBar{2}{1,1,2,2},'jb,ikac->ijkabc');    % (jk)(bc)

    DC2 =  einsum_kg(HBar{1}{1,2},L2,'kc,ijab->ijkabc') ... % (1)
         - einsum_kg(HBar{1}{1,2},L2,'ka,ijcb->ijkabc') ... % (ac)
         - einsum_kg(HBar{1}{1,2},L2,'kb,ijac->ijkabc') ... % (bc)
         - einsum_kg(HBar{1}{1,2},L2,'ic,kjab->ijkabc') ... % (ik)
         - einsum_kg(HBar{1}{1,2},L2,'jc,ikab->ijkabc') ... % (jk)
         + einsum_kg(HBar{1}{1,2},L2,'ia,kjcb->ijkabc') ... % (ik)(ac)
         + einsum_kg(HBar{1}{1,2},L2,'ib,kjac->ijkabc') ... % (ik)(bc)
         + einsum_kg(HBar{1}{1,2},L2,'ja,ikcb->ijkabc') ... % (jk)(ac)
         + einsum_kg(HBar{1}{1,2},L2,'jb,ikac->ijkabc');    % (jk)(bc)


    D1 =   einsum_kg(HBar{2}{1,2,2,2},L2,'ieab,jkec->ijkabc') ... % (1)
         - einsum_kg(HBar{2}{1,2,2,2},L2,'jeab,ikec->ijkabc') ... % (ij)
         - einsum_kg(HBar{2}{1,2,2,2},L2,'keab,jiec->ijkabc') ... % (ik)
         - einsum_kg(HBar{2}{1,2,2,2},L2,'iecb,jkea->ijkabc') ... % (ac)
         - einsum_kg(HBar{2}{1,2,2,2},L2,'ieac,jkeb->ijkabc') ... % (bc)
         + einsum_kg(HBar{2}{1,2,2,2},L2,'jecb,ikea->ijkabc') ... % (ij)(ac)
         + einsum_kg(HBar{2}{1,2,2,2},L2,'jeac,ikeb->ijkabc') ... % (ij)(bc)
         + einsum_kg(HBar{2}{1,2,2,2},L2,'kecb,jiea->ijkabc') ... % (ik)(ac)
         + einsum_kg(HBar{2}{1,2,2,2},L2,'keac,jieb->ijkabc');    % (ik)(bc)

    D2 =   einsum_kg(HBar{2}{1,1,2,1},L2,'ijam,mkbc->ijkabc') ... % (1)
         - einsum_kg(HBar{2}{1,1,2,1},L2,'kjam,mibc->ijkabc') ... % (ik)
         - einsum_kg(HBar{2}{1,1,2,1},L2,'ikam,mjbc->ijkabc') ... % (jk)
         - einsum_kg(HBar{2}{1,1,2,1},L2,'ijbm,mkac->ijkabc') ... % (ab)
         - einsum_kg(HBar{2}{1,1,2,1},L2,'ijcm,mkba->ijkabc') ... % (ac)
         + einsum_kg(HBar{2}{1,1,2,1},L2,'kjbm,miac->ijkabc') ... % (ik)(ab)
         + einsum_kg(HBar{2}{1,1,2,1},L2,'kjcm,miba->ijkabc') ... % (ik)(ac)
         + einsum_kg(HBar{2}{1,1,2,1},L2,'ikbm,mjac->ijkabc') ... % (jk)(ab)
         + einsum_kg(HBar{2}{1,1,2,1},L2,'ikcm,mjba->ijkabc');    % (jk)(ac)

    L3 = DC1 + DC2 + D1 - D2;

    fprintf(' finished in %4.2f s\n',toc)

    % permute from ijkabc to abcijk
    %L3 = permute(L3,[4,5,6,1,2,3]);

end

