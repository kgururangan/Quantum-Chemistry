function [Beta] = build_partial_MM23(i,j,k,t1,t2,h_amij,I_bcek)

    D1 =   einsum_kg(I_bcek(:,:,:,k),t2(:,:,i,j),'bce,ae->abc');
    D2 =   einsum_kg(h_amij(:,:,i,j),t2(:,:,:,k),'am,bcm->abc');

    Beta = (D1 - D2);
  
end
