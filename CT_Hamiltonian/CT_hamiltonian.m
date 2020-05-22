function [H] = CT_hamiltonian(elec_states, par_CT, par_vib)

        te = par_CT.te;
        th = par_CT.th;
        Jct = par_CT.Jct;
        num_elec_states = size(elec_states,2);
    
        n_quanta = par_vib.n_quanta;
        omega_g = par_vib.omega_g;
        Dpv_pos = par_vib.Dpv_pos;
        Dpv_neg = par_vib.Dpv_neg;
        n_modes = length(omega_g);
        num_vib_states = n_quanta^n_modes;

        % single site vibrational hamiltonian
        Ann = diag(sqrt([1:n_quanta-1]),1);
        Cre = Ann';
        Hvib_pos = zeros(num_vib_states);
        Hvib_neg = zeros(num_vib_states);
        for i = 1:n_modes
            H0temp2 = Dpv_pos(i)*omega_g(i)*(Cre+Ann);
            H0temp3 = Dpv_neg(i)*omega_g(i)*(Cre+Ann);
            A2 = 1; A3 = 1;
            for j = 1:n_modes
                if j ~= i
                    A2 = kron(eye(n_quanta),A2);
                    A3 = kron(eye(n_quanta),A3);
                else
                    A2 = kron(H0temp2,A2);
                    A3 = kron(H0temp3,A3);
                end
            end
            Hvib_pos = Hvib_pos + A2;
            Hvib_neg = Hvib_neg + A3;
        end
        

        H = zeros(num_elec_states*num_vib_states);

        for i = 1:num_elec_states
                for j = i:num_elec_states

                    A = elec_states(:,i); B = elec_states(:,j);

                    % excitation number conservation
                    if sum(A) ~= sum(B) 
                        continue
                    else

                        idx = find(A ~= B);

                        % 1-body interaction is +/- diagonal
                        if length(idx) == 2 % single excitation
                            if abs(idx(1) - idx(2)) == 1
                                kernel = zeros(num_elec_states);
                                kernel(i,j) = 1; 
                                zpq = (conj(te) - th)*eye(num_vib_states);
                                zqp = (te - conj(th))*eye(num_vib_states);
                                H = H + kron(kernel,zpq) + kron(kernel',zqp);
                            end
                        end

                        % 2-body interaction is diagonal
                        if i == j
                            
                            iocc = find(A == 1);
                            kernel = zeros(num_elec_states); kernel(i,i) = 1;
                            for k1 = 1:length(iocc)
                                for k2 = 1:length(iocc)
                                    if abs(iocc(k1) - iocc(k2)) <= par_CT.ct_range
                                        H = H + kron(kernel,-Jct*eye(num_vib_states)+Hvib_pos+Hvib_neg);
                                    end
                                end
                                H = H + kron(kernel,-Jct*eye(num_vib_states)+Hvib_pos+Hvib_neg);
                            end
                            
                        end
                        
                    end
                end
        end

end
