function [sigma] = HC_ucc_matvec(C,sys)

    sigma = zeros(size(C));
    
%     posv1 = [1:1+sys.Nocc_alpha*sys.Nvir_alpha];
%     posv2 = [posv1(end):posv1(end)+sys.Nocc_beta*sys.Nvir_beta];

    c1a = reshape(C(sys.posv{1}),sys.size{1});
    c1b = reshape(C(sys.posv{2}),sys.size{2});

    [X1A] = build_HC_1A(c1a,c1b,sys);
    [X1B] = build_HC_1B(c1a,c1b,sys);

    sigma(sys.posv{1}) = X1A(:);
    sigma(sys.posv{2}) = X1B(:);
    
end