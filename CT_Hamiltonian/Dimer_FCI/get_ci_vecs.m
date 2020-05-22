function [ci_vec] = get_ci_vecs(num_chromophore,spin_cons)

    % minimum basis set
    n_elec = 2*num_chromophore;
    norb = 2*n_elec;

    if spin_cons == 1
        
        alpha = 1:2:norb;
        beta = 2:2:norb;
        
        ci_vec = zeros(nchoosek(n_elec,length(alpha))*nchoosek(n_elec,length(beta)), n_elec);

        ci_vec_alpha = nchoosek(alpha, n_elec/2);
        ci_vec_beta = nchoosek(beta, n_elec/2);

        ct = 1;
        for i = 1:size(ci_vec_alpha,1)
            for j = 1:size(ci_vec_beta,1)
                ci_vec(ct,:) = sort([ci_vec_alpha(i,:),ci_vec_beta(j,:)],'ascend');
                ct = ct+1;
            end
        end
        
    else
        
        ci_vec = nchoosek(1:norb, n_elec);
        for i = 1:size(ci_vec,1)
            ci_vec(i,:) = sort(ci_vec(i,:),'ascend');
        end
        
    end
   
    


end

