function [sigma] = HR_ucc_matvec(R,HBar_t,cc_t,sys)

    sigma = zeros(size(R));

    r1a = reshape(R(sys.posv{1}),sys.size{1});
    r1b = reshape(R(sys.posv{2}),sys.size{2});
    r2a = reshape(R(sys.posv{3}),sys.size{3});
    r2b = reshape(R(sys.posv{4}),sys.size{4});
    r2c = reshape(R(sys.posv{5}),sys.size{5});

    [X1A] = build_HR_1A(r1a,r1b,r2a,r2b,r2c,HBar_t,cc_t,sys);
    [X1B] = build_HR_1B(r1a,r1b,r2a,r2b,r2c,HBar_t,cc_t,sys);
    [X2A] = build_HR_2A(r1a,r1b,r2a,r2b,r2c,HBar_t,cc_t,sys);
    [X2B] = build_HR_2B(r1a,r1b,r2a,r2b,r2c,HBar_t,cc_t,sys);
    [X2C] = build_HR_2C(r1a,r1b,r2a,r2b,r2c,HBar_t,cc_t,sys);

    sigma(sys.posv{1}) = X1A(:);
    sigma(sys.posv{2}) = X1B(:);
    sigma(sys.posv{3}) = X2A(:);
    sigma(sys.posv{4}) = X2B(:);
    sigma(sys.posv{5}) = X2C(:);

    
end