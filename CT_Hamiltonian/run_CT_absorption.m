function [I, omega, MatEG, omega_mat, Ve, Vg, Ee, Eg, He, Hg] = run_CT_absorption(par_CT, par_vib, par_spec)

    flag_gs_vib = 0; % no vibrational structure on gs. simply |00...0> (L sites)

    fprintf('\n')

    % build hamiltonians
    fprintf('   Building ground and CT Hamiltonians...\n')
    tic
    elec_states = get_states(par_CT.num_sites);
    [Hg] = build_gs_hamiltonian(par_CT, par_vib);
    %[He] = build_CT_Hamiltonian(elec_states, par_CT, par_vib);
    [Hf] = frenkel_hamiltonian(elec_states, par_CT, par_vib);
    %Hct = 0;
    [Hct] = CT_hamiltonian(elec_states, par_CT, par_vib);
    %Hf = 0;
    %Hct = Hct + par_CT.omega_ge*eye(size(Hct,1));
    He = Hct + Hf;
    [Ve,De] = eig(He); [Vg,Dg] = eig(Hg); 
    Ee = diag(De); Eg = diag(Dg);
    fprintf('   Hamiltonians built and diagonalized in %4.1f s\n',toc)

    %figure(1)
    subplot(221)
    scatter(1:length(Ee),Ee-par_CT.omega_ge,'MarkerFaceColor',[0,0,1])
    title('Eigenenergies')
    ylabel('E_e - \Delta_{0-0} / cm^{-1}')
    xlabel('State #')
    set(gca,'FontSize',13,'Linewidth',2,'Box','off')

    % Absorption spectrum
    [I, Id, omega, MatEG, omega_mat] = absorption_sum_over_states(Ee,Eg,Ve,Vg,par_spec);
    
    subplot(222)
    imagesc(MatEG)
    title('Transition Matrix')
    xlabel('Ground State #')
    ylabel('Excited State #')
    colormap(bluewhitered)
    colorbar
    set(gca,'FontSize',13,'Linewidth',2,'Box','off')

    %
    subplot(2,2,[3:4])
    %idx = find(Id > 0);
    plot(omega-par_CT.omega_ge,I,'Linewidth',2)
    hold on
    plot(omega-par_CT.omega_ge,Id,'r-','Linewidth',1.1)
    hold off
    title('Absorption Spectrum')
    xlabel('\omega - \Delta_{0-0} / cm^{-1}')
    ylabel('Intensity / arb.')
%     ll = legend(sprintf('Jct/Jcoul = %4.2f',par_CT.Jct/par_CT.Jcoul));
%     set(ll,'FontSize',13,'Location','NorthEast')
    set(gca,'FontSize',13,'Linewidth',2,'Box','off')

    %axis([-inf,inf,0, max(real(I))/5])
    
    
    
    fprintf('\n')
end

