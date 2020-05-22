function [I, Id, omega, MatEG, omega_mat] = absorption_sum_over_states(Ee, Eg, Ve, Vg, par_spec)

    Nw = par_spec.Nw;
    bw_Elec = par_spec.bw_Elec;
    ptol = par_spec.ptol;
    kbT = par_spec.kbT;
    gamma = par_spec.gamma;
    cent_Elec = par_spec.cent_Elec;
    popmin = par_spec.popmin;
    
    dom = bw_Elec/Nw;
    omega = (-dom*(Nw-1):2*dom:dom*(Nw-1)) + cent_Elec;

    dim_exc = size(Ve,1); dim_grd = size(Vg,1);
    
    % Thermal Density Matrix
    Dsort = diag(sort(Eg));
    Z = trace(expm(-Dsort/kbT)); 
    rho_eq = diag(expm(-Dsort/kbT)/Z);
    ind_alpha = rho_eq >= popmin;
    x0_g = rho_eq(ind_alpha);
    dim_alpha = length(x0_g);
    Ealph = Eg(ind_alpha);
    %x0_g = ones(1,dim_alpha);


    fprintf('   Calculating transition matrix...\n')
    tic
    MatEG = zeros(dim_exc, dim_alpha);
    for i = 1:dim_exc
        for j = 1:dim_alpha
            %MatEG(i,j) = conj(Ve(:,i))'*Vg(:,j); % transition matrix
            MatEG(i,j) = sum(Ve(:,i).*Vg(:,j));
        end
    end
    %MatEG = ones(size(MatEG));
    MatGE = MatEG';
    fprintf('   Transition matrix calculated in %4.1f s\n',toc)


    ind = find( abs(MatEG) >= ptol);
    [ie, ig] = ind2sub(size(MatEG),ind);
    fprintf('   %d transitions above threshold\n',length(ind))

    fprintf('   Calculating absorption spectrum...\n')
    tic
    omega_mat = zeros(length(ind),5);
    I = zeros(1,Nw); Id = zeros(1,Nw);
    for i = 1:length(ind)
        
        A = ie(i); B = ig(i);
        
        omega_mat(i,:) = [Ee(A)-Ealph(B), Ee(A), Ealph(B), MatEG(A,B)*MatGE(B,A), x0_g(B)];
        
        
        I = I + (x0_g(B)*MatEG(A,B)*MatGE(B,A)*...
                lineshape(omega,Ee(A)-Ealph(B),gamma,'gaussian'));
            
%         Id = Id + x0_g(B)*MatEG(A,B)*MatGE(B,A)*...
%                 lineshape(omega,Ee(A)-Ealph(B),1e-6,'delta');

        [~,idx] = min(abs(omega-(Ee(A)-Ealph(B))));
        Id(idx) = I(idx);
        
    end
    Id = (Id-min(I))./(max(I)-min(I));
    I = (I - min(I))./(max(I)-min(I));
    fprintf('   Absorption spectrum created in %4.1f s\n',toc)



end

