function [] = ccsdt_act_VT3_writer(subsect,spincase)

% This code writes the fully implementable (e.g. copy-paste) einsum codes
% that construct the (V*T3)_C 2-body intermediates of the B spin case
% needed for active space CCSDt codes.

    switch spincase
        
        case {'a', 'A'}
            
            switch subsect
                
                case 'vooo'
                    %  (amij)
                    m = 'm~';
                    for a = ['a','A']
                        for i = ['i','I']
                            for j  = ['j','J']
                                fprintf('\n%++++++++++++Vt3A(%s%s%s%s)++++++++++++\n',a,m,i,j)
                                out_proj = {a,m,i,j}; 
                                out_proj2 = out_proj; out_proj2{2} = 'm';
                                [d1,d2] = write_Vt3A(out_proj,'vooo');
                                write_term = @(x,y) write_einsum(x,cell2mat(out_proj2),y);
                                write_term(d1,'d1')
                                write_term(d2,'d2')
                                %[subname,sl] = get_subname(out_proj,spincase);
                                %termname = ['Vt3a.',subname];
                                fprintf('\nVt3A.%s = d1 + d2;\n',get_subname(out_proj))
                                %if flag_check
                                %    fprintf("\nfprintf('Error in %s = %4.10f\n',get_error(Vt3a.%s,%s(%s,%s,%s,%s))
                            end
                        end
                    end
                    
                case 'vvov'
                    
                    % (abie)
                    e = 'e~';
                    for a = ['a','A']
                        for b = ['b','B']
                            for i  = ['i','I']
                                fprintf('\n%++++++++++++Vt3A(%s%s%s%s)++++++++++++\n',a,b,i,e)
                                out_proj = {a,b,i,e}; 
                                out_proj2 = out_proj; out_proj2{4} = 'e';
                                [d1,d2] = write_Vt3A(out_proj,'vvov');
                                write_term = @(x,y) write_einsum(x,cell2mat(out_proj2),y);
                                write_term(d1,'d1')
                                write_term(d2,'d2')
                                fprintf('\nVt3A.%s = -d1 - d2;\n',get_subname(out_proj))
                            end
                        end
                    end
                    
                otherwise
                    disp('VT3A intermediate not supported yet')
            end
                    
        
        case {'b', 'B'}

            switch subsect

                case 'vooo'

                    % works! (amij)
                    m = 'm~';
                    for a = ['a','A']
                        for i = ['i','I']
                            for j  = ['j','J']
                                fprintf('\n%++++++++++++Vt3B(%s%s%s%s)++++++++++++\n',a,m,i,j)
                                out_proj = {a,m,i,j}; 
                                out_proj2 = out_proj; out_proj2{2} = 'm';
                                [d1,d2] = write_Vt3B(out_proj,'vooo');
                                write_term = @(x,y) write_einsum(x,cell2mat(out_proj2),y);
                                write_term(d1,'d1')
                                write_term(d2,'d2')
                                fprintf('\nVt3B.%s = d1 + d2;\n',get_subname(out_proj))
                            end
                        end
                    end

                case 'ovoo'

                    % works! (mbij)
                    m = 'm~';
                    for b = ['b','B']
                        for i = ['i','I']
                            for j  = ['j','J']
                                fprintf('\n%++++++++++++Vt3B(%s%s%s%s)++++++++++++\n',m,b,i,j)
                                out_proj = {m,b,i,j}; 
                                out_proj2 = out_proj; out_proj2{1} = 'm';
                                [d1,d2] = write_Vt3B(out_proj,'ovoo');
                                write_term = @(x,y) write_einsum(x,cell2mat(out_proj2),y);
                                write_term(d1,'d1')
                                write_term(d2,'d2')
                                fprintf('\nVt3B.%s = d1 + d2;\n',get_subname(out_proj))
                            end
                        end
                    end

                case 'vvvo'

                    % works! (abej)
                    e = 'e~';
                    for a = ['a','A']
                        for b = ['b','B']
                            for j  = ['j','J']
                                fprintf('\n%++++++++++++Vt3B(%s%s%s%s)++++++++++++\n',a,b,e,j)
                                out_proj = {a,b,e,j}; 
                                out_proj2 = out_proj; out_proj2{3} = 'e';
                                [d1,d2] = write_Vt3B(out_proj,'vvvo');
                                write_term = @(x,y) write_einsum(x,cell2mat(out_proj2),y);
                                write_term(d1,'d1')
                                write_term(d2,'d2')
                                fprintf('\nVt3B.%s = -d1 - d2;\n',get_subname(out_proj))
                            end
                        end
                    end

                case 'vvov'

                    % works! (abie)
                    e = 'e~';
                    for a = ['a','A']
                        for b = ['b','B']
                            for i  = ['i','I']
                                fprintf('\n%++++++++++++Vt3B(%s%s%s%s)++++++++++++\n',a,b,i,e)
                                out_proj = {a,b,i,e}; 
                                out_proj2 = out_proj; out_proj2{4} = 'e';
                                [d1,d2] = write_Vt3B(out_proj,'vvov');
                                write_term = @(x,y) write_einsum(x,cell2mat(out_proj2),y);
                                write_term(d1,'d1')
                                write_term(d2,'d2')
                                fprintf('\nVt3B.%s = -d1 - d2;\n',get_subname(out_proj))
                            end
                        end
                    end

                otherwise
                    disp('VT3B intermediate not supported yet')
            end
            
    end

end



function [d1,d2] = write_Vt3B(out_proj,type)

    d1 = {};
    d2 = {};
    
    switch type
        
        case 'vooo'

                % vooo (amij)
                % spin-integrated:
                % VTB.vooo = einsum_kg(sys.vB_oovv,t3b,'nmfe,afeinj->amij')...
                %            +0.5*einsum_kg(sys.vC_oovv,t3c,'nmfe,afeinj->amij');
                % n,f,e = act/unact, m = all occupied
                a = out_proj{1}; m = out_proj{2}; i = out_proj{3}; j = out_proj{4};

                % v2B * t3b
                ct = 1;
                for n = ['n','N']
                    for f = ['f','F']
                        for e = ['e','E']
                            coef = '';
                            q1 = 'v2B';
                            q2 = 't3b';
                            arr1 = {n,m,f,e};
                            arr2 = [a,f,e,i,n,j];
                            [new_arr,sign] = fix_t3_indices(arr2,q2(end));
                            term1 = [sign,coef,q1,'(',cell2mat(arr1),')'];
                            term2 = [q2,'(',new_arr,')'];
                            d1{ct} = [term1,',',term2];
                            ct = ct + 1;
                        end
                    end
                end

                % v2C * t3c
                ct = 1; 
                for n = ['n','N']
                    for f = ['f','F']
                        for e = ['e','E']

                            if f == 'f' && e == 'E' % skip redundant case
                                continue
                            end

                            val1 = get_act_idx(e); val2 = get_act_idx(f);
                            if val1 == val2
                                coef = '0.5';
                            else
                                coef = '';
                            end

                            q1 = 'v2C';
                            q2 = 't3c';
                            arr1 = {n,m,f,e};
                            arr2 = [a,f,e,i,n,j];
                            [new_arr,sign] = fix_t3_indices(arr2,q2(end));
                            term1 = [sign,coef,q1,'(',cell2mat(arr1),')'];
                            term2 = [q2,'(',new_arr,')'];
                            d2{ct} = [term1,',',term2];
                            ct = ct + 1;

                        end
                    end
                end
                
        case 'ovoo'
            
                % ovoo (mbij)
                % spin-integrated:
                % VTB.ovoo = 0.5*einsum_kg(sys.vA_oovv,t3b,'mnef,efbinj->mbij')...
                %               +einsum_kg(sys.vB_oovv,t3c,'mnef,efbinj->mbij');
             
                % n,f,e = act/unact, m = all occupied
                m = out_proj{1}; b = out_proj{2}; i = out_proj{3}; j = out_proj{4};

                % v2B * t3c
                ct = 1;
                for n = ['n','N']
                    for f = ['f','F']
                        for e = ['e','E']
                            coef = '';
                            q1 = 'v2B';
                            q2 = 't3c';
                            arr1 = {m,n,e,f};
                            arr2 = [e,f,b,i,n,j];
                            [new_arr,sign] = fix_t3_indices(arr2,q2(end));
                            term1 = [sign,coef,q1,'(',cell2mat(arr1),')'];
                            term2 = [q2,'(',new_arr,')'];
                            d1{ct} = [term1,',',term2];
                            ct = ct + 1;
                        end
                    end
                end

                % v2A * t3b
                ct = 1; 
                for n = ['n','N']
                    for f = ['f','F']
                        for e = ['e','E']

                            if f == 'f' && e == 'E' % skip redundant case
                                continue
                            end

                            val1 = get_act_idx(e); val2 = get_act_idx(f);
                            if val1 == val2
                                coef = '0.5';
                            else
                                coef = '';
                            end

                            q1 = 'v2A';
                            q2 = 't3b';
                            arr1 = {m,n,e,f};
                            arr2 = [e,f,b,i,n,j];
                            [new_arr,sign] = fix_t3_indices(arr2,q2(end));
                            term1 = [sign,coef,q1,'(',cell2mat(arr1),')'];
                            term2 = [q2,'(',new_arr,')'];
                            d2{ct} = [term1,',',term2];
                            ct = ct + 1;

                        end
                    end
                end
                
        case 'vvvo'
            
                % vvvo (abej)
                % spin-integrated:
                % VTB.vvvo = -0.5*einsum_kg(sys.vA_oovv,t3b,'mnef,afbmnj->abej')...
                %              -einsum_kg(sys.vB_oovv,t3c,'mnef,afbmnj->abej');
             
                % m,n,f = act/unact, e = all unoccupied
                a = out_proj{1}; b = out_proj{2}; e = out_proj{3}; j = out_proj{4};

                % v2B * t3c
                ct = 1;
                for m = ['m','M']
                    for n = ['n','N']
                        for f = ['f','F']
                            coef = '';
                            q1 = 'v2B';
                            q2 = 't3c';
                            arr1 = {m,n,e,f};
                            arr2 = [a,f,b,m,n,j];
                            [new_arr,sign] = fix_t3_indices(arr2,q2(end));
                            term1 = [sign,coef,q1,'(',cell2mat(arr1),')'];
                            term2 = [q2,'(',new_arr,')'];
                            d1{ct} = [term1,',',term2];
                            ct = ct + 1;
                        end
                    end
                end

                % v2A * t3b
                ct = 1; 
                for m = ['m','M']
                    for n = ['n','N']
                        for f = ['f','F']

                            if m == 'm' && n == 'N' % skip redundant case
                                continue
                            end

                            val1 = get_act_idx(m); val2 = get_act_idx(n);
                            if val1 == val2
                                coef = '0.5';
                            else
                                coef = '';
                            end

                            q1 = 'v2A';
                            q2 = 't3b';
                            arr1 = {m,n,e,f};
                            arr2 = [a,f,b,m,n,j];
                            [new_arr,sign] = fix_t3_indices(arr2,q2(end));
                            term1 = [sign,coef,q1,'(',cell2mat(arr1),')'];
                            term2 = [q2,'(',new_arr,')'];
                            d2{ct} = [term1,',',term2];
                            ct = ct + 1;

                        end
                    end
                end
                
        case 'vvov'
            
                % vvov (abie)
                % spin-integrated:
                % VTB.vvov = -einsum_kg(sys.vB_oovv,t3b,'nmfe,afbinm->abie')...
                %            -0.5*einsum_kg(sys.vC_oovv,t3c,'nmfe,afbinm->abie');
             
                % m,n,f = act/unact, e = all unoccupied
                a = out_proj{1}; b = out_proj{2}; i = out_proj{3}; e = out_proj{4};

                % v2B * t3b
                ct = 1;
                for m = ['m','M']
                    for n = ['n','N']
                        for f = ['f','F']
                            coef = '';
                            q1 = 'v2B';
                            q2 = 't3b';
                            arr1 = {n,m,f,e};
                            arr2 = [a,f,b,i,n,m];
                            [new_arr,sign] = fix_t3_indices(arr2,q2(end));
                            term1 = [sign,coef,q1,'(',cell2mat(arr1),')'];
                            term2 = [q2,'(',new_arr,')'];
                            d1{ct} = [term1,',',term2];
                            ct = ct + 1;
                        end
                    end
                end

                % v2C * t3c
                ct = 1; 
                for m = ['m','M']
                    for n = ['n','N']
                        for f = ['f','F']

                            if m == 'm' && n == 'N' % skip redundant case
                                continue
                            end

                            val1 = get_act_idx(m); val2 = get_act_idx(n);
                            if val1 == val2
                                coef = '0.5';
                            else
                                coef = '';
                            end

                            q1 = 'v2C';
                            q2 = 't3c';
                            arr1 = {n,m,f,e};
                            arr2 = [a,f,b,i,n,m];
                            [new_arr,sign] = fix_t3_indices(arr2,q2(end));
                            term1 = [sign,coef,q1,'(',cell2mat(arr1),')'];
                            term2 = [q2,'(',new_arr,')'];
                            d2{ct} = [term1,',',term2];
                            ct = ct + 1;

                        end
                    end
                end
        
            
        otherwise
            fprintf('intermediate %s not supported yet\n',type)
    end

end

function [d1,d2] = write_Vt3A(out_proj,type)

    d1 = {};
    d2 = {};
    
    switch type
        
        case 'vooo'

                % vooo (amij)
                % spin-integrated:
                % VTA.vooo = 0.5*einsum_kg(sys.vA_oovv,t3a,'mnef,aefijn->amij')...
                %              +einsum_kg(sys.vB_oovv,t3b,'mnef,aefijn->amij');
                % n,f,e = act/unact, m = all occupied
                a = out_proj{1}; m = out_proj{2}; i = out_proj{3}; j = out_proj{4};

                % v2B * t3b
                ct = 1;
                for n = ['n','N']
                    for f = ['f','F']
                        for e = ['e','E']
                            coef = '';
                            q1 = 'v2B';
                            q2 = 't3b';
                            arr1 = {m,n,e,f};
                            arr2 = [a,e,f,i,j,n];
                            [new_arr,sign] = fix_t3_indices(arr2,q2(end));
                            term1 = [sign,coef,q1,'(',cell2mat(arr1),')'];
                            term2 = [q2,'(',new_arr,')'];
                            d1{ct} = [term1,',',term2];
                            ct = ct + 1;
                        end
                    end
                end

                % v2A * t3a
                ct = 1; 
                for n = ['n','N']
                    for f = ['f','F']
                        for e = ['e','E']

                            if f == 'f' && e == 'E' % skip redundant case
                                continue
                            end

                            val1 = get_act_idx(e); val2 = get_act_idx(f);
                            if val1 == val2
                                coef = '0.5';
                            else
                                coef = '';
                            end

                            q1 = 'v2A';
                            q2 = 't3a';
                            arr1 = {m,n,e,f};
                            arr2 = [a,e,f,i,j,n];
                            [new_arr,sign] = fix_t3_indices(arr2,q2(end));
                            term1 = [sign,coef,q1,'(',cell2mat(arr1),')'];
                            term2 = [q2,'(',new_arr,')'];
                            d2{ct} = [term1,',',term2];
                            ct = ct + 1;

                        end
                    end
                end
                
        case 'vvov'
            
                % vvov (abie)
                % spin-integrated:
                % VTA.vvov = -0.5*einsum_kg(sys.vA_oovv,t3a,'mnef,abfimn->abie')...
                %              -einsum_kg(sys.vB_oovv,t3b,'mnef,abfimn->abie');
                % m,n,f = act/unact, e = all unoccupied
                a = out_proj{1}; b = out_proj{2}; i = out_proj{3}; e = out_proj{4};

                % v2B * t3b
                ct = 1;
                for m = ['m','M']
                    for n = ['n','N']
                        for f = ['f','F']
                            coef = '';
                            q1 = 'v2B';
                            q2 = 't3b';
                            arr1 = {m,n,e,f};
                            arr2 = [a,b,f,i,m,n];
                            [new_arr,sign] = fix_t3_indices(arr2,q2(end));
                            term1 = [sign,coef,q1,'(',cell2mat(arr1),')'];
                            term2 = [q2,'(',new_arr,')'];
                            d1{ct} = [term1,',',term2];
                            ct = ct + 1;
                        end
                    end
                end

                % v2A * t3a
                ct = 1; 
                for m = ['m','M']
                    for n = ['n','N']
                        for f = ['f','F']

                            if m == 'm' && n == 'N' % skip redundant case
                                continue
                            end

                            val1 = get_act_idx(m); val2 = get_act_idx(n);
                            if val1 == val2
                                coef = '0.5';
                            else
                                coef = '';
                            end

                            q1 = 'v2A';
                            q2 = 't3a';
                            arr1 = {m,n,e,f};
                            arr2 = [a,b,f,i,m,n];
                            [new_arr,sign] = fix_t3_indices(arr2,q2(end));
                            term1 = [sign,coef,q1,'(',cell2mat(arr1),')'];
                            term2 = [q2,'(',new_arr,')'];
                            d2{ct} = [term1,',',term2];
                            ct = ct + 1;

                        end
                    end
                end
        
            
        otherwise
            fprintf('intermediate %s not supported yet\n',type)
    end

end

function [new_arr,sign_char] = fix_t3_indices(arr,spin)

    sign = 1;

    switch spin
        
        case {'A','a'}
            part = arr(1:3); hole = arr(4:6);
            [perm_part] = fix_particle_idx(part);
            [perm_hole] = fix_hole_idx(hole);
            new_arr = arr([perm_part,perm_hole]);
            sign = sign * permutationParity(perm_part) * permutationParity(perm_hole);
            
        case {'B','b'}
            part = arr(1:2); hole = arr(4:5);
            [perm_part] = fix_particle_idx(part);
            [perm_hole] = fix_hole_idx(hole);
            new_arr = arr([perm_part,3,perm_hole,6]);
            sign = sign * permutationParity(perm_part) * permutationParity(perm_hole); 
            
        case {'c','C'}
            part = arr(2:3); hole = arr(5:6);
            [perm_part] = fix_particle_idx(part);
            [perm_hole] = fix_hole_idx(hole);
            new_arr = arr([1,perm_part+1,4,perm_hole+1]); % !!!
            sign = sign * permutationParity(perm_part) * permutationParity(perm_hole);

        case {'D','d'}
            part = arr(1:3); hole = arr(4:6);
            [perm_part] = fix_particle_idx(part);
            [perm_hole] = fix_hole_idx(hole);
            new_arr = arr([perm_part,perm_hole]);
            sign = sign * permutationParity(perm_part) * permutationParity(perm_hole);
            
        otherwise 
            disp('Spin type not allowed for T3!')
            
    end
    
    if sign == 1
        sign_char = '';
    else
        sign_char = '-';
    end
            

end

function [perm] = fix_particle_idx(idx)

    act_idx = zeros(1,length(idx));
    for i = 1:length(idx)
        val = get_act_idx(idx(i));
        act_idx(i) = val;
    end
    
    [~,perm] = sort(act_idx,'descend');

end

function [perm] = fix_hole_idx(idx)

    act_idx = zeros(1,length(idx));
    for i = 1:length(idx)
        val = get_act_idx(idx(i));
        act_idx(i) = val;
    end
    
    [~,perm] = sort(act_idx,'ascend');
    perm = perm + 3; %!!!
 
end

function [val] = get_character(char)

    val_act = get_act_idx(char);
    val_ph = get_ph_idx(char);
    
    if contains(char,'~') && val_ph == 1 % all occupied
        val = 5;
    elseif contains(char,'~') && val_ph == 0 % all unoccupied
        val = 6;
    else
    
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

function [name] = get_subname(idx_cell)

    name = '';
    slices = cell(1,length(idx_cell));
    for i = 1:length(idx_cell)
        
        val = get_character(idx_cell{i});
        %ph = get_ph_idx(idx_cell{i});
        switch val
            case 1
                add_char = 'P';
                %name = [name,'P'];
            case 2
                add_char = 'H';
                %name = [name,'H'];
            case 3
                add_char = 'p';
                %name = [name,'p'];
            case 4
                add_char = 'h';
                %name = [name,'h'];
            case 5
                add_char = 'v';
                %name = [name,'v'];
            case 6
                add_char = 'o';
                %name = [name,'o'];
            otherwise
                fprintf('index character %s not recognized!\n',idx_cell{i})
        end
        
        name = [name,add_char];
        
%         switch spincase
%             
%             case {'a','A'}
%                 slices{i} = [add_char,'A'];
%             case {'b','B'}
%                 if mod(i,2) == 1
%                     slices{i} = [add_char,'A'];
%                 else
%                     slices{i} = [add_char,'B'];
%                 end
%             case {'c','C'}
%                 if mod(i,2) == 1
%                     slices{i} = [add_char,'B'];
%                 else
%                     slices{i} = [add_char,'A'];
%                 end
%             case {'d','D'}
%                 slices{i} = [add_char,'B'];
%         end

    end

end


