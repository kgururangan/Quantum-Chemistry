function [dVVmat_cell] = get_dERImat(basis_cell,grad_coords)

    Nat = size(grad_coords,1);

    Norb = length(basis_cell);
    
    dVVmat_cell = cell(Nat,1);
    
    for M = 1:Nat

        dVVmat = zeros(Norb,Norb,Norb,Norb,3);
        
        RX = grad_coords(M,:);
    
        fprintf('Full derivative ERI evaluation...\n')

        tic
        % loop through permutationally unique 2-body integrals
        fprintf('ERI integrals (ij|kl)\ni = ')
        for i = 0:Norb-1
            fprintf('%d ',i+1);
            for j = 0:i
                ij = i*(i+1)/2+j;
                for k = 0:Norb-1
                    for l = 0:k
                        kl = idivide(int32(k*(k+1)),int32(2),'floor') + l;
                        if ij >= kl

                            I = i+1; J = j+1; K = k+1; L = l+1;
                            val = grad_eriABCD(basis_cell{I},basis_cell{J},basis_cell{K},basis_cell{L},RX);

                            % populate ERI tensor in chemist notation
                            dVVmat(I,J,K,L,:) = val;
                            dVVmat(J,I,K,L,:) = val;
                            dVVmat(I,J,L,K,:) = val;
                            dVVmat(K,L,I,J,:) = val;
                            dVVmat(L,K,I,J,:) = val;
                            dVVmat(L,K,J,I,:) = val;
                            dVVmat(J,I,L,K,:) = val;
                            dVVmat(K,L,J,I,:) = val;
                        end
                    end
                end
            end
        end

%         for i = 1:Norb
%             for j = 1:Norb
%                 for k = 1:Norb
%                     for l = 1:Norb
%                         dVVmat(i,j,k,l,:) = grad_eriABCD(basis_cell{i},basis_cell{j},basis_cell{k},basis_cell{l},RX);
%                     end
%                 end
%             end
%         end
        fprintf('\nFull ERI tensor evaluation took %4.2f seconds\n',toc)
        
        dVVmat_cell{M} = dVVmat;
        
    end

end