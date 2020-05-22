function [err,VVmat_fft,VVmat] = test_eri_fft(ORBS,Rat,L,h)

    Norb = length(ORBS);

    [~, ~, ~, VVmat] = spatial_integrals_v2(ORBS,Rat,ones(1,size(Rat,1)));
    
    dx = h(1); dy = h(2); dz = h(3);

    x = -L(1)/2:dx:L(1)/2; Nx = length(x);
    y = -L(2)/2:dy:L(2)/2; Ny = length(y);
    z = -L(3)/2:dz:L(3)/2; Nz = length(z);
    
    ORBS_mat = cell(1,Norb);
    for J = 1:Norb    
        
        alpha = ORBS{J}.exps;
        R = ORBS{J}.origin;
        coeff = ORBS{J}.coeff;
        
        num_gauss = length(alpha);
        
        temp = zeros(Nx,Ny,Nz);
        for i = 1:Nx
            for j = 1:Ny
                for k = 1:Nz
                    % loop over contracted primitives
                    for I = 1:num_gauss
                        temp(i,j,k) = temp(i,j,k) + primitive_norm(0,0,0,alpha(I))*...
                                                    coeff(I)*exp(-alpha(I)*((x(i)-R(1))^2+(y(j)-R(2))^2+(z(k)-R(3))^2));
                    end 
                end
            end
        end
        
        ORBS_mat{J} = temp;
    end

    
    [VVmat_fft] = eri_fft(ORBS_mat,[L,L,L],[h,h,h],0);
    
    err = zeros(Norb,Norb,Norb,Norb);
    for p = 1:Norb
        for q = 1:Norb
            for r = 1:Norb
                for s = 1:Norb
                    
                     err(p,q,r,s) = abs( (VVmat(p,q,r,s)-VVmat_fft(p,q,r,s))/VVmat(p,q,r,s) )*100;
                     
                     fprintf('Gaussian integration gives V(%d,%d,%d,%d) = %4.6f\n',p,q,r,s,VVmat(p,q,r,s))
                     fprintf('FFT integration gives V(%d,%d,%d,%d) = %4.6f\n',p,q,r,s,VVmat_fft(p,q,r,s))
                     fprintf('Error = %4.6f%% \n\n',err(p,q,r,s))
                     
                end
            end
        end
    end
    

end
