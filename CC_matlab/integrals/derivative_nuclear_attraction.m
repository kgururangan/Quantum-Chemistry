function [dxout] = derivative_nuclear_attraction(a,lmn1,A,b,lmn2,B,C,RX)

    dxout = zeros(3,1);

    if all(A == RX) && all(B == RX) && all(C == RX)
        return
    end
    
    if all(A == RX)
        
        NB = norm_factor(lmn2,b);
        NA_Lup = norm_factor(lmn1+[1,0,0],a);
        NA_Mup = norm_factor(lmn1+[0,1,0],a);
        NA_Nup = norm_factor(lmn1+[0,0,1],a);
        
        temp_x = NB*NA_Lup*sqrt((2*lmn1(1)+1)*a)*nuclear_attraction(a,[lmn1(1)+1,lmn1(2),lmn1(3)],A,b,lmn2,B,C);
        temp_y = NB*NA_Mup*sqrt((2*lmn1(2)+1)*a)*nuclear_attraction(a,[lmn1(1),lmn1(2)+1,lmn1(3)],A,b,lmn2,B,C);
        temp_z = NB*NA_Nup*sqrt((2*lmn1(3)+1)*a)*nuclear_attraction(a,[lmn1(1),lmn1(2),lmn1(3)+1],A,b,lmn2,B,C);
        
        if lmn1(1) > 0
            NA_Ldwn = norm_factor(lmn1-[1,0,0],a);
            dxout(1) = dxout(1) ...
                -NB*NA_Ldwn*2*lmn1(1)*sqrt(a/(2*lmn1(1)-1))*nuclear_attraction(a,[lmn1(1)-1,lmn1(2),lmn1(3)],A,b,lmn2,B,C);
        end
        
        if lmn1(2) > 0
            NA_Mdwn = norm_factor(lmn1-[0,1,0],a);
            dxout(2) = dxout(2) ...
                -NB*NA_Mdwn*2*lmn1(2)*sqrt(a/(2*lmn1(2)-1))*nuclear_attraction(a,[lmn1(1),lmn1(2)-1,lmn1(3)],A,b,lmn2,B,C);
        end
        
        if lmn1(3) > 0
            NA_Ndwn = norm_factor(lmn1-[0,0,1],a);
            dxout(3) = dxout(3) ...
                -NB*NA_Ndwn*2*lmn1(3)*sqrt(a/(2*lmn1(3)-1))*nuclear_attraction(a,[lmn1(1),lmn1(2),lmn1(3)-1],A,b,lmn2,B,C);
        end
        
        dxout = dxout + [temp_x; temp_y; temp_z];
    end
    
    if all(B == RX)
        
        NA = norm_factor(lmn1,a);
        NB_Lup = norm_factor(lmn2+[1,0,0],b);
        NB_Mup = norm_factor(lmn2+[0,1,0],b);
        NB_Nup = norm_factor(lmn2+[0,0,1],b);
        
        
        temp_x = NA*NB_Lup*sqrt((2*lmn2(1)+1)*b)*nuclear_attraction(a,lmn1,A,b,[lmn2(1)+1,lmn2(2),lmn2(3)],B,C);
        temp_y = NA*NB_Mup*sqrt((2*lmn2(2)+1)*b)*nuclear_attraction(a,lmn1,A,b,[lmn2(1),lmn2(2)+1,lmn2(3)],B,C);
        temp_z = NA*NB_Nup*sqrt((2*lmn2(3)+1)*b)*nuclear_attraction(a,lmn1,A,b,[lmn2(1),lmn2(2),lmn2(3)+1],B,C);
        
        if lmn2(1) > 0
            NB_Ldwn = norm_factor(lmn2-[1,0,0],b);
            dxout(1) = dxout(1) ...
                -NA*NB_Ldwn*2*lmn2(1)*sqrt(b/(2*lmn2(1)-1))*nuclear_attraction(a,lmn1,A,b,[lmn2(1)-1,lmn2(2),lmn2(3)],B,C);
        end
        
        if lmn2(2) > 0
            NB_Mdwn = norm_factor(lmn2-[0,1,0],b);
            dxout(2) = dxout(2) ...
                -NA*NB_Mdwn*2*lmn2(2)*sqrt(b/(2*lmn2(2)-1))*nuclear_attraction(a,lmn1,A,b,[lmn2(1),lmn2(2)-1,lmn2(3)],B,C);
        end
        
        if lmn2(3) > 0
            NB_Ndwn = norm_factor(lmn2-[0,0,1],b);
            dxout(3) = dxout(3) ...
                -NA*NB_Ndwn*2*lmn2(3)*sqrt(b/(2*lmn2(3)-1))*nuclear_attraction(a,lmn1,A,b,[lmn2(1),lmn2(2),lmn2(3)-1],B,C);
        end
        
        dxout = dxout + [temp_x; temp_y; temp_z];
    end
    
    if all(C == RX)
        
        NA = norm_factor(lmn1,a);
        NB = norm_factor(lmn2,b);
        
        dxout(1) = dxout(1) + NA*NB*electric_field(a,lmn1,A,b,lmn2,B,C,'x');
        dxout(2) = dxout(2) + NA*NB*electric_field(a,lmn1,A,b,lmn2,B,C,'y');
        dxout(3) = dxout(3) + NA*NB*electric_field(a,lmn1,A,b,lmn2,B,C,'z');
        
    end

end

