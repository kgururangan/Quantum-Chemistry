function [dxout] = derivative_electron_repulsion(a,lmn1,A,b,lmn2,B,c,lmn3,C,d,lmn4,D,RX)

    dxout = zeros(3,1);
    
    if all(A == RX) && all(B == RX) && all(C == RX) && all(D == RX)
        return
    end
    
    if all(A == RX)
        
        NB = norm_factor(lmn2,b);
        NC = norm_factor(lmn3,c);
        ND = norm_factor(lmn4,d);
        
        NA_Lup = norm_factor(lmn1+[1,0,0],a);
        NA_Mup = norm_factor(lmn1+[0,1,0],a);
        NA_Nup = norm_factor(lmn1+[0,0,1],a);
        
        temp_x = ND*NC*NB*NA_Lup*sqrt((2*lmn1(1)+1)*a)*electron_repulsion(a,[lmn1(1)+1,lmn1(2),lmn1(3)],A,b,lmn2,B,c,lmn3,C,d,lmn4,D);
        temp_y = ND*NC*NB*NA_Mup*sqrt((2*lmn1(2)+1)*a)*electron_repulsion(a,[lmn1(1),lmn1(2)+1,lmn1(3)],A,b,lmn2,B,c,lmn3,C,d,lmn4,D);
        temp_z = ND*NC*NB*NA_Nup*sqrt((2*lmn1(3)+1)*a)*electron_repulsion(a,[lmn1(1),lmn1(2),lmn1(3)+1],A,b,lmn2,B,c,lmn3,C,d,lmn4,D);
        
        if lmn1(1) > 0
            NA_Ldwn = norm_factor(lmn1-[1,0,0],a);
            dxout(1) = dxout(1) ...
                -ND*NC*NB*NA_Ldwn*2*lmn1(1)*sqrt(a/(2*lmn1(1)-1))*electron_repulsion(a,[lmn1(1)-1,lmn1(2),lmn1(3)],A,b,lmn2,B,c,lmn3,C,d,lmn4,D);
        end
        
        if lmn1(2) > 0
            NA_Mdwn = norm_factor(lmn1-[0,1,0],a);
            dxout(2) = dxout(2) ...
                -ND*NC*NB*NA_Mdwn*2*lmn1(2)*sqrt(a/(2*lmn1(2)-1))*electron_repulsion(a,[lmn1(1),lmn1(2)-1,lmn1(3)],A,b,lmn2,B,c,lmn3,C,d,lmn4,D);
        end
        
        if lmn1(3) > 0
            NA_Ndwn = norm_factor(lmn1-[0,0,1],a);
            dxout(3) = dxout(3) ...
                -ND*NC*NB*NA_Ndwn*2*lmn1(3)*sqrt(a/(2*lmn1(3)-1))*electron_repulsion(a,[lmn1(1),lmn1(2),lmn1(3)-1],A,b,lmn2,B,c,lmn3,C,d,lmn4,D);
        end
        
        dxout = dxout + [temp_x; temp_y; temp_z];
    end
        
    
    if all(B == RX)
        
        NA = norm_factor(lmn1,a);
        NC = norm_factor(lmn3,c);
        ND = norm_factor(lmn4,d);
        
        NB_Lup = norm_factor(lmn2+[1,0,0],b);
        NB_Mup = norm_factor(lmn2+[0,1,0],b);
        NB_Nup = norm_factor(lmn2+[0,0,1],b);
        
        
        temp_x = ND*NC*NA*NB_Lup*sqrt((2*lmn2(1)+1)*b)*electron_repulsion(a,lmn1,A,b,[lmn2(1)+1,lmn2(2),lmn2(3)],B,c,lmn3,C,d,lmn4,D);
        temp_y = ND*NC*NA*NB_Mup*sqrt((2*lmn2(2)+1)*b)*electron_repulsion(a,lmn1,A,b,[lmn2(1),lmn2(2)+1,lmn2(3)],B,c,lmn3,C,d,lmn4,D);
        temp_z = ND*NC*NA*NB_Nup*sqrt((2*lmn2(3)+1)*b)*electron_repulsion(a,lmn1,A,b,[lmn2(1),lmn2(2),lmn2(3)+1],B,c,lmn3,C,d,lmn4,D);
        
        if lmn2(1) > 0
            NB_Ldwn = norm_factor(lmn2-[1,0,0],b);
            dxout(1) = dxout(1) ...
                -ND*NC*NA*NB_Ldwn*2*lmn2(1)*sqrt(b/(2*lmn2(1)-1))*electron_repulsion(a,lmn1,A,b,[lmn2(1)-1,lmn2(2),lmn2(3)],B,c,lmn3,C,d,lmn4,D);
        end
        
        if lmn2(2) > 0
            NB_Mdwn = norm_factor(lmn2-[0,1,0],b);
            dxout(2) = dxout(2) ...
                -ND*NC*NA*NB_Mdwn*2*lmn2(2)*sqrt(b/(2*lmn2(2)-1))*electron_repulsion(a,lmn1,A,b,[lmn2(1),lmn2(2)-1,lmn2(3)],B,c,lmn3,C,d,lmn4,D);
        end
        
        if lmn2(3) > 0
            NB_Ndwn = norm_factor(lmn2-[0,0,1],b);
            dxout(3) = dxout(3) ...
                -ND*NC*NA*NB_Ndwn*2*lmn2(3)*sqrt(b/(2*lmn2(3)-1))*electron_repulsion(a,lmn1,A,b,[lmn2(1),lmn2(2),lmn2(3)-1],B,c,lmn3,C,d,lmn4,D);
        end
        
        dxout = dxout + [temp_x; temp_y; temp_z];
    end
    
    
    if all(C == RX)
        
        NA = norm_factor(lmn1,a);
        NB = norm_factor(lmn2,b);
        ND = norm_factor(lmn4,d);
        
        NC_Lup = norm_factor(lmn3+[1,0,0],c);
        NC_Mup = norm_factor(lmn3+[0,1,0],c);
        NC_Nup = norm_factor(lmn3+[0,0,1],c);
        
        
        temp_x = ND*NB*NA*NC_Lup*sqrt((2*lmn3(1)+1)*c)*electron_repulsion(a,lmn1,A,b,lmn2,B,c,[lmn3(1)+1,lmn3(2),lmn3(3)],C,d,lmn4,D);
        temp_y = ND*NB*NA*NC_Mup*sqrt((2*lmn3(2)+1)*c)*electron_repulsion(a,lmn1,A,b,lmn2,B,c,[lmn3(1),lmn3(2)+1,lmn3(3)],C,d,lmn4,D);
        temp_z = ND*NB*NA*NC_Nup*sqrt((2*lmn3(3)+1)*c)*electron_repulsion(a,lmn1,A,b,lmn2,B,c,[lmn3(1),lmn3(2),lmn3(3)+1],C,d,lmn4,D);
        
        if lmn3(1) > 0
            NC_Ldwn = norm_factor(lmn3-[1,0,0],c);
            dxout(1) = dxout(1) ...
                -ND*NB*NA*NC_Ldwn*2*lmn3(1)*sqrt(c/(2*lmn3(1)-1))*electron_repulsion(a,lmn1,A,b,lmn2,B,c,[lmn3(1)-1,lmn3(2),lmn3(3)],C,d,lmn4,D);
        end
        
        if lmn3(2) > 0
            NC_Mdwn = norm_factor(lmn3-[0,1,0],c);
            dxout(2) = dxout(2) ...
                -ND*NB*NA*NC_Mdwn*2*lmn3(2)*sqrt(c/(2*lmn3(2)-1))*electron_repulsion(a,lmn1,A,b,lmn2,B,c,[lmn3(1),lmn3(2)-1,lmn3(3)],C,d,lmn4,D);
        end
        
        if lmn3(3) > 0
            NC_Ndwn = norm_factor(lmn3-[0,0,1],c);
            dxout(3) = dxout(3) ...
                -ND*NB*NA*NC_Ndwn*2*lmn3(3)*sqrt(c/(2*lmn3(3)-1))*electron_repulsion(a,lmn1,A,b,lmn2,B,c,[lmn3(1),lmn3(2),lmn3(3)-1],C,d,lmn4,D);
        end
        
        dxout = dxout + [temp_x; temp_y; temp_z];
    end
     
    
    if all(D == RX)
        
        NA = norm_factor(lmn1,a);
        NB = norm_factor(lmn2,b);
        NC = norm_factor(lmn3,c);
        
        ND_Lup = norm_factor(lmn4+[1,0,0],d);
        ND_Mup = norm_factor(lmn4+[0,1,0],d);
        ND_Nup = norm_factor(lmn4+[0,0,1],d);
        
        
        temp_x = NC*NB*NA*ND_Lup*sqrt((2*lmn4(1)+1)*d)*electron_repulsion(a,lmn1,A,b,lmn2,B,c,lmn3,C,d,[lmn4(1)+1,lmn4(2),lmn4(3)],D);
        temp_y = NC*NB*NA*ND_Mup*sqrt((2*lmn4(2)+1)*d)*electron_repulsion(a,lmn1,A,b,lmn2,B,c,lmn3,C,d,[lmn4(1),lmn4(2)+1,lmn4(3)],D);
        temp_z = NC*NB*NA*ND_Nup*sqrt((2*lmn4(3)+1)*d)*electron_repulsion(a,lmn1,A,b,lmn2,B,c,lmn3,C,d,[lmn4(1),lmn4(2),lmn4(3)+1],D);
        
        if lmn4(1) > 0
            ND_Ldwn = norm_factor(lmn4-[1,0,0],d);
            dxout(1) = dxout(1) ...
                -NC*NB*NA*ND_Ldwn*2*lmn4(1)*sqrt(d/(2*lmn4(1)-1))*electron_repulsion(a,lmn1,A,b,lmn2,B,c,lmn3,C,d,[lmn4(1)-1,lmn4(2),lmn4(3)],D);
        end
        
        if lmn4(2) > 0
            ND_Mdwn = norm_factor(lmn4-[0,1,0],d);
            dxout(2) = dxout(2) ...
                -NC*NB*NA*ND_Mdwn*2*lmn4(2)*sqrt(d/(2*lmn4(2)-1))*electron_repulsion(a,lmn1,A,b,lmn2,B,c,lmn3,C,d,[lmn4(1),lmn4(2)-1,lmn4(3)],D);
        end
        
        if lmn4(3) > 0
            ND_Ndwn = norm_factor(lmn4-[0,0,1],d);
            dxout(3) = dxout(3) ...
                -NC*NB*NA*ND_Ndwn*2*lmn4(3)*sqrt(d/(2*lmn4(3)-1))*electron_repulsion(a,lmn1,A,b,lmn2,B,c,lmn3,C,d,[lmn4(1),lmn4(2),lmn4(3)-1],D);
        end
        
        dxout = dxout + [temp_x; temp_y; temp_z];
    end

end

