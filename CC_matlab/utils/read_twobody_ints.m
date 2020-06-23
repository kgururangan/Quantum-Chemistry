function [e2int, Vnuc] = read_twobody_ints(twobody_path,Norb)

    fid = fopen(twobody_path,'r');
    
    e2int = zeros(Norb,Norb,Norb,Norb);
    
    f_line = fgetl(fid);
    while ischar(f_line)
        
        % split line of twobody file
        f_line = split(strip(f_line));
        
        % extract indices from chemist notation
        p = str2double(f_line{1});
        r = str2double(f_line{2});
        q = str2double(f_line{3});
        s = str2double(f_line{4});
        
        % print p value for line
        %fprintf('p = %d\n',p);
        
        % twobody value (pr|qs)
        val = str2double(f_line{5});
        
        if p + q + r + s ~= 0
            % store values in e2int matrix in physics notation (pq|rs)
            e2int(p,q,r,s) = val;
        else
            % record nuclear-nuclear attraction
            Vnuc = val;
        end
        
        % get next line
        f_line = fgetl(fid);
    end
    
    fclose(fid);

end