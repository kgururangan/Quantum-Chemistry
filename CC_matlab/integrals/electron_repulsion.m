function [val] = electron_repulsion(a,lmn1,A,b,lmn2,B,c,lmn3,C,d,lmn4,D)

    l1 = lmn1(1); m1 = lmn1(2); n1 = lmn1(3);
    l2 = lmn2(1); m2 = lmn2(2); n2 = lmn2(3);
    l3 = lmn3(1); m3 = lmn3(2); n3 = lmn3(3);
    l4 = lmn4(1); m4 = lmn4(2); n4 = lmn4(3);
    [p,P] = gaussian_product_center(a,A,b,B);
    [q,Q] = gaussian_product_center(c,C,d,D);
    alpha = p*q/(p+q);
    RPQ = norm(P-Q);
    val = 0.0;
    for t = 0:l1+l2
        for u = 0:m1+m2
            for v = 0:n1+n2
                for tau = 0:l3+l4
                    for nu = 0:m3+m4
                        for phi = 0:n3+n4
                            val = val + ...
                                E(l1,l2,t,A(1)-B(1),a,b) * ...
                                E(m1,m2,u,A(2)-B(2),a,b) * ...
                                E(n1,n2,v,A(3)-B(3),a,b) * ...
                                E(l3,l4,tau,C(1)-D(1),c,d) * ...
                                E(m3,m4,nu ,C(2)-D(2),c,d) * ...
                                E(n3,n4,phi,C(3)-D(3),c,d) * ...
                                power(-1,tau+nu+phi) * ...
                                R(t+tau,u+nu,v+phi,0,alpha,P(1)-Q(1),P(2)-Q(2),P(3)-Q(3),RPQ);
                        end
                    end
                end
            end
        end
    end
    val = val*2*power(pi,2.5)/(p*q*sqrt(p+q));
end

