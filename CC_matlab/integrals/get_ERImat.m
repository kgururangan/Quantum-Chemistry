function [VVmat] = get_ERImat(varargin)

    if length(varargin) > 1
        basis_cell = varargin{1};
        pre_tol = varargin{2};
    else
        basis_cell = varargin{1};
        pre_tol = 1e-10;
    end
    
    if pre_tol == 0.0
        flag_full = true;
    else
        flag_full = false;
    end

    Norb = length(basis_cell);

    VVmat = zeros(Norb,Norb,Norb,Norb);
    
    if ~flag_full
        
        fprintf('ERI evaluation with pre-screening...\n')
        %%% integral pre-screening
        % In chemist notation, the ERI tensor V(a,b,c,d) can be shown to be a positive
        % definite matrix if we take the Norb^2 x Norb^2 matrix V(ab,cd) where
        % a,b are orbitals of electron 1 and c,d are orbitals of electron 2.
        % Using the Cauchy-Schwarz inequality, 
        % |V(ab,cd)| <= sqrt(V(ab,ab))*sqrt(V(cd,cd)) therefore prescreening is achieved
        % by first calculating the diagonal elements V(ab,ab)

        % DIAGONAL ERI EVALUATION
        tic_diag = tic;

        for i = 0:Norb-1
            for j = 0:i
                ij = i*(i+1)/2+j;
                for k = i
                    for l = j
                        kl = idivide(int32(k*(k+1)),int32(2),'floor') + l;
                        if ij >= kl

                            I = i+1; J = j+1; K = k+1; L = l+1;
                            val = eriABCD(basis_cell{I},basis_cell{J},basis_cell{K},basis_cell{L});

                            % populate ERI tensor in chemist notation
                            VVmat(I,J,K,L) = val;
                            VVmat(J,I,K,L) = val;
                            VVmat(I,J,L,K) = val;
                            VVmat(K,L,I,J) = val;
                            VVmat(L,K,I,J) = val;
                            VVmat(L,K,J,I) = val;
                            VVmat(J,I,L,K) = val;
                            VVmat(K,L,J,I) = val;

                            %fprintf('\nVVmat(%d,%d,%d,%d) = %4.10f',I,J,K,L,val);
                        end
                    end
                end
            end
        end
        fprintf('Diagonal ERI tensor evaluation took %4.2f seconds\n',toc(tic_diag))

        % SCREENED 2-BODY ERI EVALUATION
        tic_pre = tic;
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

                            %if i,j,k,l -> diagonal element
                            
                            %else
                                pre = sqrt(VVmat(I,J,I,J))*sqrt(VVmat(K,L,K,L));

                                if pre > pre_tol 
                                    val = eriABCD(basis_cell{I},basis_cell{J},basis_cell{K},basis_cell{L});
                                else
                                    val = 0.0;
                                end
                             %end


                            % populate ERI tensor in chemist notation
                            VVmat(I,J,K,L) = val;
                            VVmat(J,I,K,L) = val;
                            VVmat(I,J,L,K) = val;
                            VVmat(K,L,I,J) = val;
                            VVmat(L,K,I,J) = val;
                            VVmat(L,K,J,I) = val;
                            VVmat(J,I,L,K) = val;
                            VVmat(K,L,J,I) = val;

                            %fprintf('\nVVmat(%d,%d,%d,%d) = %4.10f',I,J,K,L,val);
                        end
                    end
                end
            end
        end
        fprintf('\nScreened ERI tensor evaluation took %4.2f seconds\n',toc(tic_pre))
    
    else
        
        fprintf('Full ERI evaluation...\n')
        % FULL 2-BODY ERI EVALUATION
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
                            val = eriABCD(basis_cell{I},basis_cell{J},basis_cell{K},basis_cell{L});

                            % populate ERI tensor in chemist notation
                            VVmat(I,J,K,L) = val;
                            VVmat(J,I,K,L) = val;
                            VVmat(I,J,L,K) = val;
                            VVmat(K,L,I,J) = val;
                            VVmat(L,K,I,J) = val;
                            VVmat(L,K,J,I) = val;
                            VVmat(J,I,L,K) = val;
                            VVmat(K,L,J,I) = val;

                            %fprintf('\nVVmat(%d,%d,%d,%d) = %4.10f',I,J,K,L,val);
                        end
                    end
                end
            end
        end
        fprintf('\nFull ERI tensor evaluation took %4.2f seconds\n',toc)
        
    end

end
