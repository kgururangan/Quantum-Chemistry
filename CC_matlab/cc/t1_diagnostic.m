function [val] = t1_diagnostic(cc_t,sys)

    t1 = cat(1,cc_t.t1a(:),cc_t.t1b(:));
    N = sys.Nelec;
    
    val = norm(t1)/sqrt(N);

end

