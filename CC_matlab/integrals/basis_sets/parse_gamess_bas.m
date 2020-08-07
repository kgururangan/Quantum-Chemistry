function [basis_cell] = parse_gamess_bas(bas_file)

    fid = fopen(bas_file,'r');
    
    tline = fgetl(fid);
    flag = false;
    while ~flag
        if strcmp(tline,'$DATA')
            flag = true;
        else
            tline = fgetl(fid);
        end
    end

    tline = fgetl(fid);
    while ischar(tline) 
        xbool = is_shell(tline);
        if ~xbool
            if length(split(tline)) > 1
                L = split(tline);
                prim_exps(ct) = str2double(L{2});
                prim_coef(ct) = str2double(L{3});
                ct = ct + 1;
            end
        else
            ct = 1;
            [~, num_prim, prim_shell] = parse_shell_line(tline);
            prim_cell = cell(1,length(prim_shell));
            prim_exps = zeros(1,num_prim);
            prim_coef = zeros(1,num_prim);
        end
            
        tline = fgetl(fid);
    end
    
    fclose(fid);
    
    
end

function [xbool] = is_shell(line)

    line = split(line); shell_type = line{1};
    
    if strcmp(shell_type,'S') || strcmp(shell_type,'P') || strcmp(shell_type,'D') || strcmp(shell_type,'F') || strcmp(shell_type,'L') 
        
        xbool = true; 
        
    else
        xbool = false;
    end
    
end


function [is_shell, num_prim, prim_shell] = parse_shell_line(line)

    line = split(line); shell_type = line{1}; num_prim = NaN; prim_shell = cell(1,1);
    
    if strcmp(shell_type,'S') || strcmp(shell_type,'P') || strcmp(shell_type,'D') || strcmp(shell_type,'F') || strcmp(shell_type,'L')
        
        is_shell = true; num_prim = str2double(line{2});
        
        if strcmp(shell_type,'S')
            prim_shell = {[0, 0, 0]};
        elseif strcmp(shell_type,'L')
            prim_shell = {[0,0,0],[1,0,0],[0,1,0],[0,0,1]};
        elseif strcmp(shell_type,'P')
            prim_shell = {[1,0,0],[0,1,0],[0,0,1]};
        elseif strcmp(shell_type,'D')
            prim_shell = {[2,0,0],[0,2,0],[0,0,2],[1,1,0],[1,0,1],[0,1,1]};
        elseif strcmp(shell_type,'F')
            prim_shell = {[3,0,0],[0,3,0],[0,0,3],[2,1,0],[2,0,1],[1,2,0],[1,0,2],[0,1,2],[0,2,1],[1,1,1]};
        else
            fprintf('\nShell type %s not supported!',shell_type)
        end
        
    else
        is_shell = false;
    end
end

