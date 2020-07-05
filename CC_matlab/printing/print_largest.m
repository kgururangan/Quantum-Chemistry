function [] = print_largest(t,nprint)

    ndim = length(size(t)); tflat = t(:);
    [~,idx] = sort(abs(tflat),'descend');
    
    fprintf('Largest %d amplitudes are:\n',nprint)
    
    for p = 1:nprint
        if ndim == 4
            % NOT GOOD ORDER FOR LAMBDA VECTORS!
            [a,b,i,j] = ind2sub(size(t),idx(p));
            fprintf('t(%d,%d,%d,%d) = %4.8f\n',a,b,i,j,tflat(idx(p)));
        else
            [a,i] = ind2sub(size(t),idx(p));
            fprintf('t(%d,%d) = %4.8f\n',a,i,tflat(idx(p)));
        end
    end
    
end