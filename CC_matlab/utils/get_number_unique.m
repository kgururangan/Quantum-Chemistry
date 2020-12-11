function [num_unique] = get_number_unique(spintype,sys)

    switch spintype
        
        case '1A'
            num_unique = sys.Nocc_alpha * sys.Nvir_alpha;
        case '1B'
            num_unique = sys.Nocc_beta * sys.Nvir_beta;
        case '2A'
            num_unique = 0;
            for a = 1:sys.Nvir_alpha
                for b = a+1:sys.Nvir_alpha
                    for i = 1:sys.Nocc_alpha
                        for j = i+1:sys.Nocc_alpha
                            num_unique = num_unique + 1;
                        end
                    end
                end
            end
        case '2B'
            num_unique = 0;
            for a = 1:sys.Nvir_alpha
                for b = 1:sys.Nvir_beta
                    for i = 1:sys.Nocc_alpha
                        for j = 1:sys.Nocc_beta
                            num_unique = num_unique + 1;
                        end
                    end
                end
            end
        case '2C'
            num_unique = 0;
            for a = 1:sys.Nvir_beta
                for b = a+1:sys.Nvir_beta
                    for i = 1:sys.Nocc_beta
                        for j = i+1:sys.Nocc_beta
                            num_unique = num_unique + 1;
                        end
                    end
                end
            end
        case '3A'
            num_unique = 0;
            for a = 1:sys.Nvir_alpha
                for b = a+1:sys.Nvir_alpha
                    for c = b+1:sys.Nvir_alpha
                        for i = 1:sys.Nocc_alpha
                            for j = i+1:sys.Nocc_alpha
                                for k = j+1:sys.Nocc_alpha
                                    num_unique = num_unique + 1;
                                end
                            end
                        end
                    end
                end
            end
        case '3B'
            num_unique = 0;
            for a = 1:sys.Nvir_alpha
                for b = a+1:sys.Nvir_alpha
                    for c = 1:sys.Nvir_beta
                        for i = 1:sys.Nocc_alpha
                            for j = i+1:sys.Nocc_alpha
                                for k = 1:sys.Nocc_beta
                                    num_unique = num_unique + 1;
                                end
                            end
                        end
                    end
                end
            end
        case '3C'
            num_unique = 0;
            for a = 1:sys.Nvir_alpha
                for b = 1:sys.Nvir_beta
                    for c = b+1:sys.Nvir_beta
                        for i = 1:sys.Nocc_alpha
                            for j = 1:sys.Nocc_beta
                                for k = j+1:sys.Nocc_beta
                                    num_unique = num_unique + 1;
                                end
                            end
                        end
                    end
                end
            end
        case '3D'
            num_unique = 0;
            for a = 1:sys.Nvir_beta
                for b = a+1:sys.Nvir_beta
                    for c = b+1:sys.Nvir_beta
                        for i = 1:sys.Nocc_beta
                            for j = i+1:sys.Nocc_beta
                                for k = j+1:sys.Nocc_beta
                                    num_unique = num_unique + 1;
                                end
                            end
                        end
                    end
                end
            end
        case '4A'
            num_unique = 0;
            for a = 1:sys.Nvir_alpha
                for b = a+1:sys.Nvir_alpha
                    for c = b+1:sys.Nvir_alpha
                        for d = c+1:sys.Nvir_alpha
                            for i = 1:sys.Nocc_alpha
                                for j = i+1:sys.Nocc_alpha
                                    for k = j+1:sys.Nocc_alpha
                                        for l = k+1:sys.Nocc_alpha
                                            num_unique = num_unique + 1;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        case '4B'
            num_unique = 0;
            for a = 1:sys.Nvir_alpha
                for b = a+1:sys.Nvir_alpha
                    for c = b+1:sys.Nvir_alpha
                        for d = 1:sys.Nvir_beta
                            for i = 1:sys.Nocc_alpha
                                for j = i+1:sys.Nocc_alpha
                                    for k = j+1:sys.Nocc_alpha
                                        for l = 1:sys.Nocc_beta
                                            num_unique = num_unique + 1;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        case '4C'
            num_unique = 0;
            for a = 1:sys.Nvir_alpha
                for b = a+1:sys.Nvir_alpha
                    for c = 1:sys.Nvir_beta
                        for d = c+1:sys.Nvir_beta
                            for i = 1:sys.Nocc_alpha
                                for j = i+1:sys.Nocc_alpha
                                    for k = 1:sys.Nocc_beta
                                        for l = k+1:sys.Nocc_beta
                                            num_unique = num_unique + 1;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        case '4D'
            num_unique = 0;
            for a = 1:sys.Nvir_alpha
                for b = 1:sys.Nvir_beta
                    for c = b+1:sys.Nvir_beta
                        for d = c+1:sys.Nvir_beta
                            for i = 1:sys.Nocc_alpha
                                for j = 1:sys.Nocc_beta
                                    for k = j+1:sys.Nocc_beta
                                        for l = k+1:sys.Nocc_beta
                                            num_unique = num_unique + 1;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        case '4E'
            num_unique = 0;
            for a = 1:sys.Nvir_beta
                for b = a+1:sys.Nvir_beta
                    for c = b+1:sys.Nvir_beta
                        for d = c+1:sys.Nvir_beta
                            for i = 1:sys.Nocc_beta
                                for j = i+1:sys.Nocc_beta
                                    for k = j+1:sys.Nocc_beta
                                        for l = k+1:sys.Nocc_beta
                                            num_unique = num_unique + 1;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
        otherwise
            fprintf('Spin case %s not supported yet.\n',spintype)
    end
    
end

