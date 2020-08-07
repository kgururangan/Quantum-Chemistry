function [Beta] = build_partial_MM23(i,j,k,t1,t2,h_amij,I_bcek)

    D1 =   einsum_kg(I_bcek(:,:,:,k),t2(:,:,i,j),'bce,ae->abc') ...
         - einsum_kg(I_bcek(:,:,:,k),t2(:,:,i,j),'ace,be->abc') ...
         - einsum_kg(I_bcek(:,:,:,k),t2(:,:,i,j),'bae,ce->abc');

    D2 =   einsum_kg(h_amij(:,:,i,j),t2(:,:,:,k),'am,bcm->abc') ...
         - einsum_kg(h_amij(:,:,i,j),t2(:,:,:,k),'bm,acm->abc') ...
         - einsum_kg(h_amij(:,:,i,j),t2(:,:,:,k),'cm,bam->abc');

    Beta = D1 - D2;
  
end

