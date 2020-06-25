function [Gamma] = build_partial_L3(i,j,k,L1,L2,HBar)

    DC1 = einsum_kg(L1(k,:),HBar{2}{1,1,2,2}(i,j,:,:),'c,ab->abc');

    DC2 = einsum_kg(L2(i,j,:,:),HBar{1}{1,2}(k,:),'ab,c->abc');

    D1 = einsum_kg(L2(:,k,:,:),HBar{2}{1,1,2,1}(i,j,:,:),'mbc,am->abc');

    D2 = einsum_kg(L2(i,j,:,:),HBar{2}{2,1,2,2}(:,k,:,:),'ae,ebc->abc');

    Gamma = 0.5*(DC1 + DC2 - D1 + D2);
end
