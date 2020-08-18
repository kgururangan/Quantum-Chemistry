function [val] = electric_field(a,lmn1,A,b,lmn2,B,C,direction)

    l1 = lmn1(1); m1 = lmn1(2); n1 = lmn1(3);
    l2 = lmn2(1); m2 = lmn2(2); n2 = lmn2(3);
    [p,P] = gaussian_product_center(a,A,b,B);
    RPC = norm(P-C);
    val = 0.0;
    
    switch direction
        
        case 'x'
            
            for t = 0:l1+l2
                for u = 0:m1+m2
                    for v = 0:n1+n2
                        val = val + ...
                            E(l1,l2,t,A(1)-B(1),a,b)*E(m1,m2,u,A(2)-B(2),a,b)*E(n1,n2,v,A(3)-B(3),a,b)...
                            * R(t+1,u,v,0,p,P(1)-C(1),P(2)-C(2),P(3)-C(3),RPC);
                    end
                end
            end
            
        case 'y'
            
           for t = 0:l1+l2
                for u = 0:m1+m2
                    for v = 0:n1+n2
                        val = val + ...
                            E(l1,l2,t,A(1)-B(1),a,b)*E(m1,m2,u,A(2)-B(2),a,b)*E(n1,n2,v,A(3)-B(3),a,b)...
                            * R(t,u+1,v,0,p,P(1)-C(1),P(2)-C(2),P(3)-C(3),RPC);
                    end
                end
           end
            
        case 'z'
            
           for t = 0:l1+l2
                for u = 0:m1+m2
                    for v = 0:n1+n2
                        val = val + ...
                            E(l1,l2,t,A(1)-B(1),a,b)*E(m1,m2,u,A(2)-B(2),a,b)*E(n1,n2,v,A(3)-B(3),a,b)...
                            * R(t,u,v+1,0,p,P(1)-C(1),P(2)-C(2),P(3)-C(3),RPC);
                    end
                end
           end
           
        otherwise
            
            disp('Enter a valid Cartesian direction for the electric field integral!')
            
    end
    
    val = val * (-2*pi/p);

end
