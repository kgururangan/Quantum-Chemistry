function [H] = build_CT_Hamiltonian(elec_states, par_CT, par_vib)

    te = par_CT.te;
    th = par_CT.th;
    Jct = par_CT.Jct;
    Jcoul = par_CT.Jcoul;
    num_sites = par_CT.num_sites;
    omega_ge = par_CT.omega_ge;
    num_elec_states = size(elec_states,2);
    
    if ~isempty(par_vib)
    
            n_quanta = par_vib.n_quanta;
            omega_g = par_vib.omega_g;
            Dpv0 = par_vib.Dpv0;
            Dpv_pos = par_vib.Dpv_pos;
            Dpv_neg = par_vib.Dpv_neg;

            n_modes = length(omega_g);

            num_vib_states = n_quanta^n_modes;

            % single site vibrational hamiltonian
            Ann = diag(sqrt([1:n_quanta-1]),1);
            Cre = Ann';
            Num = Cre*Ann;

            Hvib0 = zeros(num_vib_states);
            Hvib_pos = zeros(num_vib_states);
            Hvib_neg = zeros(num_vib_states);
            for i = 1:n_modes

                H0temp1 = omega_g(i)*Num + Dpv0(i)*omega_g(i)*(Cre+Ann);
                H0temp2 = Dpv_pos(i)*omega_g(i)*(Cre+Ann);
                H0temp3 = Dpv_neg(i)*omega_g(i)*(Cre+Ann);

                A1 = 1; A2 = 1; A3 = 1;
                for j = 1:n_modes
                    if j ~= i
                        A1 = kron(eye(n_quanta),A1);
                        A2 = kron(eye(n_quanta),A2);
                        A3 = kron(eye(n_quanta),A3);
                    else
                        A1 = kron(H0temp1,A1);
                        A2 = kron(H0temp2,A2);
                        A3 = kron(H0temp3,A3);
                    end
                end
                Hvib0 = Hvib0 + A1;
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
                                % phase due to re-ordering for slater's rules
                                % here..?
                                kernel = zeros(num_elec_states);
                                kernel(i,j) = 1; 
                                zpq = (conj(te) - th)*eye(num_vib_states);
                                zqp = (te - conj(th))*eye(num_vib_states);
                                H = H + kron(kernel,zpq) + kron(kernel',zqp);
                            end
                        end

                        % 2-body interaction is diagonal
                        if i == j
                            iocc = find(elec_states(:,i) == 1);
                            for k1 = 1:length(iocc)
                                for k2 = 1:length(iocc)
                                    if abs(iocc(k1) - iocc(k2)) <= num_sites

                                        kernel = zeros(num_elec_states); kernel(i,i) = 1;
                                        
                                        H = H + kron(kernel,Jct*eye(num_vib_states)+Hvib_pos+Hvib_neg);

                                        if abs(iocc(k1) - iocc(k2)) == 1
                                            H = H + kron(kernel,Jcoul*eye(num_vib_states));
                                        end
                                        
                                        if abs(iocc(k1)-iocc(k2)) == 0
                                            H = H + kron(kernel,Hvib0);
                                        end
                                            

                                    end
                                end
                            end
                        end

                    end


                end
            end

            H = H + omega_ge*eye(size(H,1));



    else


        H = zeros(num_elec_states);

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
                            % phase due to re-ordering for slater's rules
                            % here..?
                            H(i,j) = (conj(te) - th);
                            H(j,i) = (te - conj(th));
                        end
                    end

                    % 2-body interaction is diagonal
                    if i == j
                        iocc = find(elec_states(:,i) == 1);
                        for k1 = 1:length(iocc)
                            for k2 = 1:length(iocc)
                                if abs(iocc(k1) - iocc(k2)) <= num_sites

                                    H(i,i) = H(i,i) + Jct;

                                    if abs(iocc(k1) - iocc(k2)) == 1
                                        H(i,i) = H(i,i) - Jcoul;
                                    end

                                end
                            end
                        end
                    end

                end


            end
        end

        H = H + omega_ge*eye(size(H,1));

    end

    
end

     