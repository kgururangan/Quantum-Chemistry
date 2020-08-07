function [Rnew] = update_R(R,HBar_t,sys,omega,shift)

    Nocc_a = sys.Nocc_alpha; Nocc_b = sys.Nocc_beta; 
    Nunocc_a = sys.Nvir_alpha; Nunocc_b = sys.Nvir_beta;
    
    H1A = HBar_t.H1A; H1B = HBar_t.H1B;

    r1a = reshape(R(sys.posv{1}),sys.size{1});
    r1b = reshape(R(sys.posv{2}),sys.size{2});
    r2a = reshape(R(sys.posv{3}),sys.size{3});
    r2b = reshape(R(sys.posv{4}),sys.size{4});
    r2c = reshape(R(sys.posv{5}),sys.size{5});

    r1a_new = zeros(Nunocc_a, Nocc_a);
    for i = 1:Nocc_a
        for a = 1:Nunocc_a
            denom = H1A.vv(a,a) - H1A.oo(i,i);
            r1a_new(a,i) = r1a(a,i)/(omega - denom + shift); 
        end
    end
    
    r1b_new = zeros(Nunocc_b,Nocc_b);
    for i = 1:Nocc_b
        for a = 1:Nunocc_b
            denom = H1B.vv(a,a) - H1B.oo(i,i);
            r1b_new(a,i) = r1b(a,i)/(omega - denom + shift); 
        end
    end
    
    r2a_new = zeros(Nunocc_a, Nunocc_a, Nocc_a, Nocc_a);
    for i = 1:Nocc_a
        for a = 1:Nunocc_a
            for j = i+1:Nocc_a
                for b = a+1:Nunocc_a
                    
                    denom = H1A.vv(a,a) - H1A.oo(i,i) + H1A.vv(b,b) - H1A.oo(j,j);
                      
                    r2a_new(a,b,i,j) = r2a(a,b,i,j)/(omega - denom + shift);                   
                    r2a_new(b,a,i,j) = -r2a_new(a,b,i,j);
                    r2a_new(a,b,j,i) = -r2a_new(a,b,i,j);
                    r2a_new(b,a,j,i) = r2a_new(a,b,i,j);
                end
            end
        end
    end
    
   r2b_new = zeros(Nunocc_a, Nunocc_b, Nocc_a, Nocc_b);
    for i = 1:Nocc_a
        for a = 1:Nunocc_a
            for j = 1:Nocc_b
                for b = 1:Nunocc_b                    
                    denom = H1A.vv(a,a) - H1A.oo(i,i) + H1B.vv(b,b) - H1B.oo(j,j);                      
                    r2b_new(a,b,i,j) = r2b(a,b,i,j)/(omega - denom + shift);                   
                end
            end
        end
    end
    
    r2c_new = zeros(Nunocc_b, Nunocc_b, Nocc_b, Nocc_b);
    for i = 1:Nocc_b
        for a = 1:Nunocc_b
            for j = i+1:Nocc_b
                for b = a+1:Nunocc_b
                    
                    denom = H1B.vv(a,a) - H1B.oo(i,i) + H1B.vv(b,b) - H1B.oo(j,j);
                      
                    r2c_new(a,b,i,j) = r2c(a,b,i,j)/(omega - denom + shift);                   
                    r2c_new(b,a,i,j) = -r2c_new(a,b,i,j);
                    r2c_new(a,b,j,i) = -r2c_new(a,b,i,j);
                    r2c_new(b,a,j,i) = r2c_new(a,b,i,j);
                end
            end
        end
    end

    Rnew = cat(1,r1a_new(:),r1b_new(:),r2a_new(:),r2b_new(:),r2c_new(:));


end