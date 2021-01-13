function [] = write_einsum(str_array,out_proj,term_label)


    if nargin < 3
        term_label = '';
    end

    for j = 1:length(str_array)
        
        flag = false;
        string = str_array{j};
        
        sign = string(1);
        if sign ~= '+' && sign ~= '-'
            sign = '+';
            flag = true;
        end
        
        for jj = 1:length(string)
            if string(jj) == 'h' || string(jj) == 'v' || string(jj) == 'f'
                break
            end
        end
        
        if jj == 1 && sign ~= '-'
            coef = '';
        elseif jj == 2 && sign == '-'
            coef = '';
        else
            if ~flag % positive without +
                coef = string(2:jj-1);
            else % sign included
                coef = string(1:jj-1);
            end
        end
        
        s = split(string(jj:end),',');

        num_terms = length(s);
        
        if num_terms == 2
            
            term1  = s{1}; term2 = s{2};
            s1 = split(term1,'('); q1 = s1{1}; c1 = s1{2}(1:end-1);
            s2 = split(term2,'('); q2 = s2{1}; c2 = s2{2}(1:end-1);
            
            C1 = {}; c1pr = '';
            ct = 1;
            for jjj = 1:length(c1)
                if isletter(c1(jjj))
                    C1{ct} = c1(jjj);
                    c1pr = [c1pr,c1(jjj)];
                    ct = ct + 1;
                else
                    C1{ct-1} = [C1{ct-1},c1(jjj)];
                end
            end
            
            C2 = {}; c2pr = '';
            ct = 1;
            for jjj = 1:length(c2)
                if isletter(c2(jjj))
                    C2{ct} = c2(jjj);
                    c2pr = [c2pr,c2(jjj)];
                    ct = ct + 1;
                else
                    C2{ct-1} = [C2{ct-1},c2(jjj)];
                end
            end
            
            switch q1
                
                case {'h1a','h1A','H1a','H1A'}
                    
                    q1 = 'H1A';
                    
                    if length(C1) > 2
                        disp('Error: number of indices does not match object!')
                    end

   
                    hbar_str = get_HBar_str(C1);
                    arr = get_slicing_arrays(C1);
                    
                    q1 = [q1,'.',hbar_str];
                    
                    % THIS PATCH PERMITS ONLY ONE ~ (full range) INDEX PER
                    % TERM!!!
                    ARR = [arr{1},'A',',',arr{2},'A'];
                    ARR2 = ARR;
                    for pp = 1:length(ARR)
                        if ARR(pp) == ':'
                            ARR2(pp+1) = [];
                        end
                    end
                            
                    q1 = [q1,'(',ARR2,')'];    
                      
                case {'f1a','f1A','F1a','F1A'}
                    
                    q1 = 'sys.fA';
                    
                    if length(C1) > 2
                        disp('Error: number of indices does not match object!')
                    end

                    hbar_str = get_HBar_str(C1);
                    arr = get_slicing_arrays(C1);
                    
                    q1 = [q1,'_',hbar_str];
                    ARR = [arr{1},'A',',',arr{2},'A'];
                    ARR2 = ARR;
                    for pp = 1:length(ARR)
                        if ARR(pp) == ':'
                            ARR2(pp+1) = [];
                        end
                    end
                            
                    q1 = [q1,'(',ARR2,')'];       
                    
                case {'h1b','h1B','H1b','H1B'}
                    
                    q1 = 'H1B';
                    
                    if length(C1) > 2
                        disp('Error: number of indices does not match object!')
                    end

                    hbar_str = get_HBar_str(C1);
                    arr = get_slicing_arrays(C1);
                    
                    q1 = [q1,'.',hbar_str];
                    ARR = [arr{1},'B',',',arr{2},'B'];
                    ARR2 = ARR;
                    for pp = 1:length(ARR)
                        if ARR(pp) == ':'
                            ARR2(pp+1) = [];
                        end
                    end
                            
                    q1 = [q1,'(',ARR2,')']; 
                      
                case {'f1b','f1B','F1b','F1B'}
                    
                    q1 = 'sys.fB';
                    
                    if length(C1) > 2
                        disp('Error: number of indices does not match object!')
                    end

                    hbar_str = get_HBar_str(C1);
                    arr = get_slicing_arrays(C1);
                    
                    q1 = [q1,'_',hbar_str];
                    ARR = [arr{1},'B',',',arr{2},'B'];
                    ARR2 = ARR;
                    for pp = 1:length(ARR)
                        if ARR(pp) == ':'
                            ARR2(pp+1) = [];
                        end
                    end
                            
                    q1 = [q1,'(',ARR2,')']; 
                    
                case {'h2a','h2A','H2a','H2A'}
                    
                    q1 = 'H2A';
                    
                    if length(C1) > 4
                        disp('Error: number of indices does not match object!')
                    end

                    hbar_str = get_HBar_str(C1);
                    arr = get_slicing_arrays(C1);
                    
                    q1 = [q1,'.',hbar_str];
                    ARR = [arr{1},'A',',',arr{2},'A',',',arr{3},'A',',',arr{4},'A'];
                    ARR2 = ARR;
                    for pp = 1:length(ARR)
                        if ARR(pp) == ':'
                            ARR2(pp+1) = [];
                        end
                    end
                            
                    q1 = [q1,'(',ARR2,')']; 
                      
                case {'v2a','v2A','V2a','V2A'}
                    
                    q1 = 'sys.vA';
                    
                    if length(C1) > 4
                        disp('Error: number of indices does not match object!')
                    end

                    hbar_str = get_HBar_str(C1);
                    arr = get_slicing_arrays(C1);
                    
                    q1 = [q1,'_',hbar_str];
                    ARR = [arr{1},'A',',',arr{2},'A',',',arr{3},'A',',',arr{4},'A'];
                    ARR2 = ARR;
                    for pp = 1:length(ARR)
                        if ARR(pp) == ':'
                            ARR2(pp+1) = [];
                        end
                    end
                            
                    q1 = [q1,'(',ARR2,')']; 
                    
                case {'h2b','h2B','H2b','H2B'}
                    
                    q1 = 'H2B';
                    
                    if length(C1) > 4
                        disp('Error: number of indices does not match object!')
                    end

                    hbar_str = get_HBar_str(C1);
                    arr = get_slicing_arrays(C1);
                    
                    q1 = [q1,'.',hbar_str];
                    ARR = [arr{1},'A',',',arr{2},'B',',',arr{3},'A',',',arr{4},'B'];
                    ARR2 = ARR;
                    for pp = 1:length(ARR)
                        if ARR(pp) == ':'
                            ARR2(pp+1) = [];
                        end
                    end
                            
                    q1 = [q1,'(',ARR2,')']; 
                      
                case {'v2b','v2B','V2b','V2B'}
                    
                    q1 = 'sys.vB';
                    
                    if length(C1) > 4
                        disp('Error: number of indices does not match object!')
                    end

                    hbar_str = get_HBar_str(C1);
                    arr = get_slicing_arrays(C1);
                    
                    q1 = [q1,'_',hbar_str];
                    
                    ARR = [arr{1},'A',',',arr{2},'B',',',arr{3},'A',',',arr{4},'B'];
                    ARR2 = ARR;
                    for pp = 1:length(ARR)
                        if ARR(pp) == ':'
                            ARR2(pp+1) = [];
                        end
                    end
                            
                    q1 = [q1,'(',ARR2,')']; 
                    
                case {'h2c','h2C','H2c','H2C'}   
                    
                    q1 = 'H2C';
                    
                    if length(C1) > 4
                        disp('Error: number of indices does not match object!')
                    end

                    hbar_str = get_HBar_str(C1);
                    arr = get_slicing_arrays(C1);
                    
                    q1 = [q1,'.',hbar_str];
                    ARR = [arr{1},'B',',',arr{2},'B',',',arr{3},'B',',',arr{4},'B'];
                    ARR2 = ARR;
                    for pp = 1:length(ARR)
                        if ARR(pp) == ':'
                            ARR2(pp+1) = [];
                        end
                    end
                            
                    q1 = [q1,'(',ARR2,')']; 
                      
                case {'v2c','v2C','V2c','V2C'}   
                    
                    q1 = 'sys.vC';
                    
                    if length(C1) > 4
                        disp('Error: number of indices does not match object!')
                    end

                    hbar_str = get_HBar_str(C1);
                    arr = get_slicing_arrays(C1);
                    
                    q1 = [q1,'_',hbar_str];
                    ARR = [arr{1},'B',',',arr{2},'B',',',arr{3},'B',',',arr{4},'B'];
                    ARR2 = ARR;
                    for pp = 1:length(ARR)
                        if ARR(pp) == ':'
                            ARR2(pp+1) = [];
                        end
                    end
                            
                    q1 = [q1,'(',ARR2,')']; 
                      
                otherwise
                    disp('Case for term 1 not supported!')
                    
            end
            
            
            switch q2
                
                case {'t3a','t3A','T3A','T3a'}
                    
                    q2 = 'T3A';
                    
                    arr = get_slicing_arrays(C2);
                    
                    q2 = [q2,'.',...
                          arr{1},arr{2},arr{3},arr{4},arr{5},arr{6}];
                          
                    
                case {'t3b','t3B','T3B','T3b'}
                    
                    q2 = 'T3B';
                    
                    arr = get_slicing_arrays(C2);
                    
                    q2 = [q2,'.',...
                          arr{1},arr{2},arr{3},arr{4},arr{5},arr{6}];
                    
                case {'t3c','t3C','T3C','T3c'}
                    
                    q2 = 'T3C';
                    
                    arr = get_slicing_arrays(C2);
                    
                    q2 = [q2,'.',...
                          arr{1},arr{2},arr{3},arr{4},arr{5},arr{6}];
                        
                case {'t3d','t3D','T3D','T3d'}
                    
                    q2 = 'T3D';
                    
                    arr = get_slicing_arrays(C2);
                    
                    q2 = [q2,'.',...
                          arr{1},arr{2},arr{3},arr{4},arr{5},arr{6}];
                    
            end 
            
            if ~isempty(coef)
            
                if j == 1
                    fprintf("\n%s = %s%s*einsum_kg(%s,%s,'%s,%s->%s')...",...
                             term_label,sign,coef,q1,q2,c1pr,c2pr,out_proj);
                elseif j < length(str_array)
                    fprintf("\n%s%s*einsum_kg(%s,%s,'%s,%s->%s')...",...
                             sign,coef,q1,q2,c1pr,c2pr,out_proj);
                else
                    fprintf("\n%s%s*einsum_kg(%s,%s,'%s,%s->%s');\n",...
                             sign,coef,q1,q2,c1pr,c2pr,out_proj);
                end
                
            else
                
                if j == 1
                    fprintf("\n%s = %seinsum_kg(%s,%s,'%s,%s->%s')...",...
                             term_label,sign,q1,q2,c1pr,c2pr,out_proj);
                elseif j < length(str_array)
                    fprintf("\n%seinsum_kg(%s,%s,'%s,%s->%s')...",...
                             sign,q1,q2,c1pr,c2pr,out_proj);
                else
                    fprintf("\n%seinsum_kg(%s,%s,'%s,%s->%s');\n",...
                             sign,q1,q2,c1pr,c2pr,out_proj);
                end
                
            end
            
            
        end
            
    end
    


end


function [arr] = get_slicing_arrays(indices_cell)

    arr = cell(1,length(indices_cell));
    for i = 1:length(indices_cell)
        val = get_character(indices_cell{i});
        if val == 0
            arr{i} = ':';
        end
        if val == 1
            arr{i} = 'P';
        end
        if val == 2
            arr{i} = 'H';
        end
        if val == 3
            arr{i} = 'p';
        end
        if val == 4
            arr{i} = 'h';
        end
    end

end


function [str] = get_HBar_str(indices_cell)

    str = '';

    for i = 1:length(indices_cell)
        val = get_ph_idx(indices_cell{i});
        if val == 1
            str = [str,'v'];
        end
        if val == 0
            str = [str,'o'];
        end
    end

end
    

function [val] = get_character(char)

    val_act = get_act_idx(char);
    val_ph = get_ph_idx(char);
    
    if contains(char,'~')% && val_ph == 1 % all occupied
        val = 0;
    else
    
%     if contains(char,'~') && val_ph == 0 % all unoccupied
%         val = 6;
%     end
    
        if val_act == 1 && val_ph == 1 % active particle
            val = 1;
        end

        if val_act == 1 && val_ph == 0 % active hole
            val = 2;
        end

        if val_act == 0 && val_ph == 1 % virtual
            val = 3;
        end

        if val_act == 0 && val_ph == 0 % core
            val = 4;
        end
        
    end


end

function [val] = get_act_idx(char)

    if length(char) > 1
        char1 = char(1);
    else
        char1 = char;
    end

    if isstrprop(char1,'upper')
        val = 1;
    else
        val = 0;
    end
end

function [val] = get_ph_idx(char)


    if length(char) > 1
        char1 = char(1);
    else
        char1 = char;
    end

    switch char1
        
        case {'m','n','i','j','k','M','N','I','J','K'}
            val = 0;
            
        case {'e','f','a','b','c','E','F','A','B','C'}
            val = 1;
            
    end
              
end
