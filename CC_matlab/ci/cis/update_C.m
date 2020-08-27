function [Cnew] = update_C(C,sys,omega,shift)

    Nocc_a = sys.Nocc_alpha; Nocc_b = sys.Nocc_beta; 
    Nunocc_a = sys.Nvir_alpha; Nunocc_b = sys.Nvir_beta;
    
    c1a = reshape(C(sys.posv{1}),sys.size{1});
    c1b = reshape(C(sys.posv{2}),sys.size{2});

    c1a_new = zeros(Nunocc_a, Nocc_a);
    for i = 1:Nocc_a
        for a = 1:Nunocc_a
            %denom = sys.fa_vv(a,a) - sys.fa_oo(i,i);
            denom = sys.fa_vv(a,a) - sys.fa_oo(i,i) + sys.vA_ovvo(i,a,a,i);
            c1a_new(a,i) = c1a(a,i)/(omega - denom + shift); 
        end
    end
    
    c1b_new = zeros(Nunocc_b,Nocc_b);
    for i = 1:Nocc_b
        for a = 1:Nunocc_b
            %denom = sys.fb_vv(a,a) - sys.fb_oo(i,i);
            denom = sys.fb_vv(a,a) - sys.fb_oo(i,i) + sys.vC_ovvo(i,a,a,i);
            c1b_new(a,i) = c1b(a,i)/(omega - denom + shift); 
        end
    end
%     
%     r2a_new = zeros(Nunocc_a, Nunocc_a, Nocc_a, Nocc_a);
%     for i = 1:Nocc_a
%         for a = 1:Nunocc_a
%             for j = i+1:Nocc_a
%                 for b = a+1:Nunocc_a
%                     
%                     denom = H1A.vv(a,a) - H1A.oo(i,i) + H1A.vv(b,b) - H1A.oo(j,j);
%                       
%                     r2a_new(a,b,i,j) = r2a(a,b,i,j)/(omega - denom + shift);                   
%                     r2a_new(b,a,i,j) = -r2a_new(a,b,i,j);
%                     r2a_new(a,b,j,i) = -r2a_new(a,b,i,j);
%                     r2a_new(b,a,j,i) = r2a_new(a,b,i,j);
%                 end
%             end
%         end
%     end
%     
%    r2b_new = zeros(Nunocc_a, Nunocc_b, Nocc_a, Nocc_b);
%     for i = 1:Nocc_a
%         for a = 1:Nunocc_a
%             for j = 1:Nocc_b
%                 for b = 1:Nunocc_b                    
%                     denom = H1A.vv(a,a) - H1A.oo(i,i) + H1B.vv(b,b) - H1B.oo(j,j);                      
%                     r2b_new(a,b,i,j) = r2b(a,b,i,j)/(omega - denom + shift);                   
%                 end
%             end
%         end
%     end
%     
%     r2c_new = zeros(Nunocc_b, Nunocc_b, Nocc_b, Nocc_b);
%     for i = 1:Nocc_b
%         for a = 1:Nunocc_b
%             for j = i+1:Nocc_b
%                 for b = a+1:Nunocc_b
%                     
%                     denom = H1B.vv(a,a) - H1B.oo(i,i) + H1B.vv(b,b) - H1B.oo(j,j);
%                       
%                     r2c_new(a,b,i,j) = r2c(a,b,i,j)/(omega - denom + shift);                   
%                     r2c_new(b,a,i,j) = -r2c_new(a,b,i,j);
%                     r2c_new(a,b,j,i) = -r2c_new(a,b,i,j);
%                     r2c_new(b,a,j,i) = r2c_new(a,b,i,j);
%                 end
%             end
%         end
%     end

    Cnew = cat(1,c1a_new(:),c1b_new(:));


end