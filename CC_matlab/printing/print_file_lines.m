function [] = print_file_lines(fid)

    tline = fgetl(fid);
    while ischar(tline)
        disp(tline)
        tline = fgetl(fid);
    end
    fclose(fid);
    
end


