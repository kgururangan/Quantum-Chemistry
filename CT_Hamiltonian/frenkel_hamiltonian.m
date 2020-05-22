function [H] = frenkel_hamiltonian(elec_states, par_CT, par_vib)

        Jcoul = par_CT.Jcoul;
        omega_ge = par_CT.omega_ge;
        num_elec_states = size(elec_states,2);
    
        n_quanta = par_vib.n_quanta;
        omega_g = par_vib.omega_g;
        Dpv0 = par_vib.Dpv0;
        n_modes = length(omega_g);
        num_vib_states = n_quanta^n_modes;

        % single site vibrational hamiltonian
        Ann = diag(sqrt([1:n_quanta-1]),1);
        Cre = Ann';
        Num = Cre*Ann;

        Hvib0 = zeros(num_vib_states);
        for i = 1:n_modes

            H0temp1 = omega_g(i)*Num + Dpv0(i)*omega_g(i)*(Cre+Ann);

            A1 = 1; 
            for j = 1:n_modes
                if j ~= i
                    A1 = kron(eye(n_quanta),A1);
                else
                    A1 = kron(H0temp1,A1);
                end
            end
            Hvib0 = Hvib0 + A1;
        end
        

        H = zeros(num_elec_states*num_vib_states);

        for i = 1:num_elec_states

                iocc = find(elec_states(:,i) == 1);
                
                kernel = zeros(num_elec_states); kernel(i,i) = 1;

                for k1 = 1:length(iocc)

                    %kernel = zeros(num_elec_states); kernel(iocc(k1),iocc(k1)) = 1;


                    for k2 = 1:length(iocc)
                        if abs(iocc(k1) - iocc(k2)) == 1
                            H = H + kron(kernel,2*Jcoul*eye(num_vib_states));
                        end
                    end

                    H = H + kron(kernel,Hvib0);

                end

        end

        H = H + omega_ge*eye(size(H,1));

end




