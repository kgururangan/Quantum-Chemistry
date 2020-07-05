function [e1int, Norb] = read_onebody_ints(varargin)

    fid = fopen(varargin{1},'r');
    
    if length(varargin) == 1
        num_lines = 0;
        f_line = fgetl(fid);
        while ischar(f_line)
            num_lines = num_lines + 1;
            f_line = fgetl(fid);
        end
        frewind(fid);

        Norb = 0.5*(-1 + sqrt(1+8*num_lines));
    else
        Norb = varargin{2};
    end
    
    e1int = zeros(Norb);
    for p = 1:Norb
        for q = 1:p
            f_line = split(strip(fgetl(fid)));
            e1int(p,q) = str2double(f_line{1});
            e1int(q,p) = str2double(f_line{1});
        end
    end

    fclose(fid);

end

