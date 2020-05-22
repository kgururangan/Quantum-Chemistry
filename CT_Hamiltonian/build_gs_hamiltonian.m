function [Hg] = build_gs_hamiltonian(par_CT, par_vib)

    if ~isempty(par_vib) 
        
        n_quanta = par_vib.n_quanta;
        omega_g = par_vib.omega_g;
        G = par_vib.anharmonic_tensor;

        n_modes = length(omega_g);

        N = n_quanta^n_modes;

        % single site vibrational hamiltonian
        Ann = diag(sqrt([1:n_quanta-1]),1);
        Cre = Ann';
        Num = Cre*Ann;


        Htemp = zeros(N,N);
        for i = 1:n_modes

            H0temp = omega_g(i)*Num;

            A = 1;
            for j = 1:n_modes
                if j ~= i
                    A = kron(eye(n_quanta),A);
                else
                    A = kron(H0temp,A);
                end
            end
            Htemp = Htemp + A;
        end


        % Anharmonic correction to ground state
        Hanh_g = zeros(N,N); % ground-state harmonic Hamiltonian

        % Hanh_g = kron(kron(Ann + Cre,Ann + Cre),Ann + Cre);
        for i = 1:n_modes
            for j = 1:n_modes
                for k = 1:n_modes
                    vec = [i,j,k];
                    ALL = 1;
                    for l = 1:n_modes
                       [Il,Jl] = find(vec == l);
                       nl = length(Jl);
                       ALL = kron((Ann + Cre)^nl,ALL);
                    end
                    Hanh_g = Hanh_g + G(i,j,k)*ALL;
                end
            end
        end

        Hanh_g = (1/6)*Hanh_g;

        Hanh_g = Hanh_g + Htemp;

        %Hg = kron(eye(2^par_CT.num_sites),Hanh_g);
        Hg = kron(Hanh_g, eye(2^par_CT.num_sites));

    else
        
        Hg = 0;
        
    end


end

