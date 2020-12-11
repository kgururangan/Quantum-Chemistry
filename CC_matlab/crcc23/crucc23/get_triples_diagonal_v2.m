function [d3a1,d3a2,d3b1,d3b2,d3c1,d3c2,d3d1,d3d2] = get_triples_diagonal_v2(cc_t,sys)

    fprintf('\nCalculating 3-body triples diagonal... ')
    tic

    Nocc_a = sys.Nocc_alpha; Nunocc_a = sys.Nvir_alpha;
    Nocc_b = sys.Nocc_beta; Nunocc_b = sys.Nvir_beta;

    t2a = cc_t.t2a; t2b = cc_t.t2b; t2c = cc_t.t2c;

    d3a1 = zeros(Nunocc_a,Nocc_a,Nocc_a);
    for i = 1:Nocc_a
        for j = 1:Nocc_a
            for a = 1:Nunocc_a
                val = 0.0;
                for e = 1:Nunocc_a
                    val = val + sys.vA_vvoo(a,e,i,j) * t2a(a,e,i,j);
                end
                d3a1(a,j,i) = val;
            end
        end
    end
    
    d3a2 = zeros(Nunocc_a,Nunocc_a,Nocc_a);
    for i = 1:Nocc_a
        for a = 1:Nunocc_a
            for b = 1:Nunocc_a
                val = 0.0;
                for m = 1:Nocc_a
                    val = val + sys.vA_vvoo(a,b,i,m) * t2a(a,b,i,m);
                end
                d3a2(b,a,i) = -val;
            end
        end
    end
    
    d3b1 = zeros(Nunocc_a,Nocc_b,Nocc_a);
    for i = 1:Nocc_a
        for j = 1:Nocc_b
            for a = 1:Nunocc_a
                val = 0.0;
                for e = 1:Nunocc_b
                    val = val + sys.vB_vvoo(a,e,i,j) * t2b(a,e,i,j);
                end
                d3b1(a,j,i) = val;
            end
        end
    end
    
    d3b2 = zeros(Nunocc_b,Nunocc_a,Nocc_a);
    for i = 1:Nocc_a
        for a = 1:Nunocc_a
            for b = 1:Nunocc_b
                val = 0.0;
                for m = 1:Nocc_b
                    val = val + sys.vB_vvoo(a,b,i,m) * t2b(a,b,i,m);
                end
                d3b2(b,a,i) = -val;
            end
        end
    end
    
    d3c1 = zeros(Nunocc_b,Nocc_b,Nocc_a);
    for i = 1:Nocc_a
        for j = 1:Nocc_b
            for b = 1:Nunocc_b
                val = 0.0;
                for e = 1:Nunocc_a
                    val = val + sys.vB_vvoo(e,b,i,j) * t2b(e,b,i,j);
                end
                d3c1(b,j,i) = val;
            end
        end
    end
    
    d3c2 = zeros(Nunocc_b,Nunocc_a,Nocc_b);
    for j = 1:Nocc_b
        for a = 1:Nunocc_a
            for b = 1:Nunocc_b
                val = 0.0;
                for m = 1:Nocc_a
                    val = val + sys.vB_vvoo(a,b,m,j) * t2b(a,b,m,j);
                end
                d3c2(b,a,j) = -val;
            end
        end
    end
    
    d3d1 = zeros(Nunocc_b,Nocc_b,Nocc_b);
    for i = 1:Nocc_b
        for j = 1:Nocc_b
            for a = 1:Nunocc_b
                val = 0.0;
                for e = 1:Nunocc_b
                    val = val + sys.vC_vvoo(a,e,i,j) * t2c(a,e,i,j);
                end
                d3d1(a,j,i) = val;
            end
        end
    end
    
    d3d2 = zeros(Nunocc_b,Nunocc_b,Nocc_b);
    for i = 1:Nocc_b
        for a = 1:Nunocc_b
            for b = 1:Nunocc_b
                val = 0.0;
                for m = 1:Nocc_b
                    val = val + sys.vC_vvoo(a,b,i,m) * t2c(a,b,i,m);
                end
                d3d2(b,a,i) = -val;
            end
        end
    end
    

    fprintf('finished in %4.4f s\n',toc);
    
end
