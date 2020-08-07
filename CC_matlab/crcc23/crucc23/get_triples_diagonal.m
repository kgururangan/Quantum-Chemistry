function [D3A_V,D3A_O,D3B_V,D3B_O,D3C_V,D3C_O,D3D_V,D3D_O] = get_triples_diagonal(cc_t,sys)

    fprintf('\nCalculating 3-body triples diagonal... ')
    tic

    Nocc_a = sys.Nocc_alpha; Nunocc_a = sys.Nvir_alpha;
    Nocc_b = sys.Nocc_beta; Nunocc_b = sys.Nvir_beta;

    t2a = cc_t.t2a; t2b = cc_t.t2b; t2c = cc_t.t2c;

    % D3A2(b,a,i) -> 
    %d3A_V = @(a,i,b) -dot(squeeze(sys.vA_oovv(i,:,b,a)),squeeze(t2a(a,b,i,:))); % works better but wrong..?
    d3A_V = @(a,i,b) -dot(squeeze(sys.vA_oovv(i,:,a,b)),squeeze(t2a(a,b,i,:)));
    % D3A1(a,j,i) -> 
    d3A_O = @(a,i,j) dot(squeeze(sys.vA_oovv(i,j,a,:)),squeeze(t2a(a,:,i,j)));
    
    % D3B2(b,a,i) -> b = c, 
    d3B_V = @(a,i,c) -dot(squeeze(sys.vB_oovv(i,:,a,c)),squeeze(t2b(a,c,i,:)));
    % D3B1(a,j,i) -> j = k
    d3B_O = @(a,i,k) dot(squeeze(sys.vB_oovv(i,k,a,:)),squeeze(t2b(a,:,i,k)));
    
    % D3C2(b,a,j) -> j = k, b = c
    d3C_V = @(a,k,c) -dot(squeeze(sys.vB_oovv(:,k,a,c)),squeeze(t2b(a,c,:,k)));
    % D3C1(b,j,i) -> j = k, b = c
    d3C_O = @(c,i,k) dot(squeeze(sys.vB_oovv(i,k,:,c)),squeeze(t2b(:,c,i,k)));
    
    % D3D2(b,a,i) -> 
    %d3D_V = @(a,i,b) -dot(squeeze(sys.vC_oovv(i,:,b,a)),squeeze(t2c(a,b,i,:)));
    d3D_V = @(a,i,b) -dot(squeeze(sys.vC_oovv(i,:,a,b)),squeeze(t2c(a,b,i,:)));
    % D3D1(a,j,i) ->
    d3D_O = @(a,i,j) dot(squeeze(sys.vC_oovv(i,j,a,:)),squeeze(t2c(a,:,i,j)));
    
    D3A_V = zeros(Nunocc_a,Nocc_a,Nunocc_a);
    D3A_O = zeros(Nunocc_a,Nocc_a,Nocc_a);
    D3B_V = zeros(Nunocc_a,Nocc_a,Nunocc_b);
    D3B_O = zeros(Nunocc_a,Nocc_a,Nocc_b);
    D3C_V = zeros(Nunocc_a,Nocc_b,Nunocc_b);
    D3C_O = zeros(Nunocc_b,Nocc_a,Nocc_b);
    D3D_V = zeros(Nunocc_b,Nocc_b,Nunocc_b);
    D3D_O = zeros(Nunocc_b,Nocc_b,Nocc_b);


    
    % A diagonal
    for a = 1:Nunocc_a
        for i = 1:Nocc_a
            for b = 1:Nunocc_a
                % D3A2(b,a,i)
                D3A_V(a,i,b) = d3A_V(a,i,b);
            end
        end
    end
    
    for a = 1:Nunocc_a
        for i = 1:Nocc_a
            for j = 1:Nocc_a
                % D3A1(a,j,i)
                D3A_O(a,i,j) = d3A_O(a,i,j);
            end
        end
    end
    
    % B diagonal
    for a = 1:Nunocc_a
        for i = 1:Nocc_a
            for c = 1:Nunocc_b
                % D3B2(b,a,i)
                D3B_V(a,i,c) = d3B_V(a,i,c);
            end
        end
    end
    
    for a = 1:Nunocc_a
        for i = 1:Nocc_a
            for k = 1:Nocc_b
                % D3B1(a,j,i)
                D3B_O(a,i,k) = d3B_O(a,i,k);
            end
        end
    end
    
   % C diagonal 
    for a = 1:Nunocc_a
        for k = 1:Nocc_b
            for c = 1:Nunocc_b
                % D3C2(b,a,j)
                D3C_V(a,k,c) = d3C_V(a,k,c);
            end
        end
    end
    
    for c = 1:Nunocc_b
        for i = 1:Nocc_a
            for k = 1:Nocc_b
                % D3C1(b,j,i)
                D3C_O(c,i,k) = d3C_O(c,i,k);
            end
        end
    end
    
    % D diagonal 
    for a = 1:Nunocc_b
        for i = 1:Nocc_b
            for b = 1:Nunocc_b
                % D3D2(b,a,i)
                D3D_V(a,i,b) = d3D_V(a,i,b);
            end
        end
    end
    
    for a = 1:Nunocc_b
        for i = 1:Nocc_b
            for j = 1:Nocc_b
                % D3D1(a,j,i)
                D3D_O(a,i,j) = d3D_O(a,i,j);
            end
        end
    end
    
    fprintf('finished in %4.4f s\n',toc);
    
end

