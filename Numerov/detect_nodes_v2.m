function [ZC,NN] = detect_nodes_v2(y,dir)
    
    if isempty(y)
        ZC = []; NN = 1e3;
        return 
    end

    NN = 0; ZC = [];
    if strcmp(dir,'left')
        for i = 1:length(y)-1
            if (y(i) > 0 && y(i+1) <0) || (y(i) < 0 && y(i+1) > 0) || (y(i) == 0)
                NN = NN + 1;
                ZC(NN) = i;
            end
        end
    else
        for i = 2:length(y)
            if (y(i) > 0 && y(i-1) <0) || (y(i) < 0 && y(i-1) > 0) || (y(i) == 0)
                NN = NN + 1;
                ZC(NN) = i;
            end
        end
    end
    
end

