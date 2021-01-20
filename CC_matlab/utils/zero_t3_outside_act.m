function [x] = zero_t3_outside_act(x_in,num_act_thresh,spin_case,sys)

    x = x_in;
    
    switch spin_case

        case {'A','a'}

                ctP = 0;
                ctQ = 0;

                Nunocc_a = sys.Nvir_alpha;
                Nocc_a = sys.Nocc_alpha;

                act_h_alpha_range = [sys.Nocc_alpha - sys.Nact_h_alpha + 1 : sys.Nocc_alpha];
                act_p_alpha_range = [1 : sys.Nact_p_alpha];

                for a = 1:Nunocc_a
                    for b = a+1:Nunocc_a
                        for c = b+1:Nunocc_a
                            for i = 1:Nocc_a
                                for j = i+1:Nocc_a
                                    for k = j+1:Nocc_a
                                        num_act_h = num_active([i,j,k],act_h_alpha_range);
                                        num_act_p = num_active([a,b,c],act_p_alpha_range);
                                        if num_act_h >=num_act_thresh && num_act_p >= num_act_thresh
                                            ctP = ctP + 1;
                                            continue
                                        else          
                                            x(a,b,c,i,j,k) = 0;
                                            x(a,b,c,k,i,j) = x(a,b,c,i,j,k);
                                            x(a,b,c,j,k,i) = x(a,b,c,i,j,k);
                                            x(a,b,c,i,k,j) = -x(a,b,c,i,j,k);
                                            x(a,b,c,j,i,k) = -x(a,b,c,i,j,k);
                                            x(a,b,c,k,j,i) = -x(a,b,c,i,j,k);

                                            % (ab)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                            x(b,a,c,i,j,k) = -x(a,b,c,i,j,k);
                                            x(b,a,c,k,i,j) = -x(a,b,c,i,j,k);
                                            x(b,a,c,j,k,i) = -x(a,b,c,i,j,k);
                                            x(b,a,c,i,k,j) = x(a,b,c,i,j,k);
                                            x(b,a,c,j,i,k) = x(a,b,c,i,j,k);
                                            x(b,a,c,k,j,i) = x(a,b,c,i,j,k);

                                            % (bc)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                            x(a,c,b,i,j,k) = -x(a,b,c,i,j,k);
                                            x(a,c,b,k,i,j) = -x(a,b,c,i,j,k);
                                            x(a,c,b,j,k,i) = -x(a,b,c,i,j,k);
                                            x(a,c,b,i,k,j) = x(a,b,c,i,j,k);
                                            x(a,c,b,j,i,k) = x(a,b,c,i,j,k);
                                            x(a,c,b,k,j,i) = x(a,b,c,i,j,k);

                                            % (ac)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                            x(c,b,a,i,j,k) = -x(a,b,c,i,j,k);
                                            x(c,b,a,k,i,j) = -x(a,b,c,i,j,k);
                                            x(c,b,a,j,k,i) = -x(a,b,c,i,j,k);
                                            x(c,b,a,i,k,j) = x(a,b,c,i,j,k);
                                            x(c,b,a,j,i,k) = x(a,b,c,i,j,k);
                                            x(c,b,a,k,j,i) = x(a,b,c,i,j,k);

                                            % (ac)(bc)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                            x(b,c,a,i,j,k) = x(a,b,c,i,j,k);
                                            x(b,c,a,k,i,j) = x(a,b,c,i,j,k);
                                            x(b,c,a,j,k,i) = x(a,b,c,i,j,k);
                                            x(b,c,a,i,k,j) = -x(a,b,c,i,j,k);
                                            x(b,c,a,j,i,k) = -x(a,b,c,i,j,k);
                                            x(b,c,a,k,j,i) = -x(a,b,c,i,j,k);

                                            % (ac)(ab)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                            x(c,a,b,i,j,k) = x(a,b,c,i,j,k);
                                            x(c,a,b,k,i,j) = x(a,b,c,i,j,k);
                                            x(c,a,b,j,k,i) = x(a,b,c,i,j,k);
                                            x(c,a,b,i,k,j) = -x(a,b,c,i,j,k);
                                            x(c,a,b,j,i,k) = -x(a,b,c,i,j,k);
                                            x(c,a,b,k,j,i) = -x(a,b,c,i,j,k);
                                            ctQ = ctQ + 1;
                                        end

                                    end
                                end
                            end
                        end
                    end
                end

            case {'B','b'}

                ctP = 0;
                ctQ = 0;

                Nunocc_a = sys.Nvir_alpha;
                Nunocc_b = sys.Nvir_beta;
                Nocc_a = sys.Nocc_alpha;
                Nocc_b = sys.Nocc_beta;

                act_h_alpha_range = [sys.Nocc_alpha - sys.Nact_h_alpha + 1 : sys.Nocc_alpha];
                act_h_beta_range = [sys.Nocc_beta - sys.Nact_h_beta + 1 : sys.Nocc_beta];
                act_p_alpha_range = [1 : sys.Nact_p_alpha];
                act_p_beta_range = [1 : sys.Nact_p_beta];

                for a = 1:Nunocc_a
                    for b = a+1:Nunocc_a
                        for c = 1:Nunocc_b
                            for i = 1:Nocc_a
                                for j = i+1:Nocc_a
                                    for k = 1:Nocc_b
                                        num_act_h = num_active([i,j],act_h_alpha_range) + num_active([k],act_h_beta_range);
                                        num_act_p = num_active([a,b],act_p_alpha_range) + num_active([c],act_p_beta_range);

                                        if num_act_h >=num_act_thresh && num_act_p >= num_act_thresh
                                            ctP = ctP + 1;
                                            continue
                                        else          
                                            x(a,b,c,i,j,k) = 0;
                                            x(b,a,c,i,j,k) = -x(a,b,c,i,j,k);
                                            x(a,b,c,j,i,k) = -x(a,b,c,i,j,k);
                                            x(b,a,c,j,i,k) = x(a,b,c,i,j,k);
                                            ctQ = ctQ + 1;
                                        end

                                        end
                                    end
                                end
                            end
                        end
                    end
                

        case {'C','c'}

                ctP = 0;
                ctQ = 0;

                Nunocc_a = sys.Nvir_alpha;
                Nunocc_b = sys.Nvir_beta;
                Nocc_a = sys.Nocc_alpha;
                Nocc_b = sys.Nocc_beta;

                act_h_alpha_range = [sys.Nocc_alpha - sys.Nact_h_alpha + 1 : sys.Nocc_alpha];
                act_h_beta_range = [sys.Nocc_beta - sys.Nact_h_beta + 1 : sys.Nocc_beta];
                act_p_alpha_range = [1 : sys.Nact_p_alpha];
                act_p_beta_range = [1 : sys.Nact_p_beta];

                for a = 1:Nunocc_a
                    for b = 1:Nunocc_b
                        for c = b+1:Nunocc_b
                            for i = 1:Nocc_a
                                for j = 1:Nocc_b
                                    for k = j+1:Nocc_b

                                        num_act_h = num_active([i],act_h_alpha_range) + num_active([j,k],act_h_beta_range);
                                        num_act_p = num_active([a],act_p_alpha_range) + num_active([b,c],act_p_beta_range);

                                        if num_act_h >=num_act_thresh && num_act_p >= num_act_thresh
                                            ctP = ctP + 1;
                                            continue
                                        else          
                                            x(a,b,c,i,j,k) = 0;
                                            x(a,c,b,i,j,k) = -x(a,b,c,i,j,k);
                                            x(a,b,c,i,k,j) = -x(a,b,c,i,j,k);
                                            x(a,c,b,i,k,j) = x(a,b,c,i,j,k);
                                            ctQ = ctQ + 1;
                                        end

                                    end
                                end
                            end
                        end
                    end
                end
            

        case {'D','d'}

                ctP = 0;
                ctQ = 0;

                Nunocc_b = sys.Nvir_beta;
                Nocc_b = sys.Nocc_beta;

                act_h_beta_range = [sys.Nocc_beta - sys.Nact_h_beta + 1 : sys.Nocc_beta];
                act_p_beta_range = [1 : sys.Nact_p_beta];


                for a = 1:Nunocc_b
                    for b = a+1:Nunocc_b
                        for c = b+1:Nunocc_b
                            for i = 1:Nocc_b
                                for j = i+1:Nocc_b
                                    for k = j+1:Nocc_b

                                        num_act_h = num_active([i,j,k],act_h_beta_range);
                                        num_act_p = num_active([a,b,c],act_p_beta_range);

                                        if num_act_h >=num_act_thresh && num_act_p >= num_act_thresh
                                            ctP = ctP + 1;
                                            continue
                                        else          
                                            x(a,b,c,i,j,k) = 0;
                                            x(a,b,c,k,i,j) = x(a,b,c,i,j,k);
                                            x(a,b,c,j,k,i) = x(a,b,c,i,j,k);
                                            x(a,b,c,i,k,j) = -x(a,b,c,i,j,k);
                                            x(a,b,c,j,i,k) = -x(a,b,c,i,j,k);
                                            x(a,b,c,k,j,i) = -x(a,b,c,i,j,k);

                                            % (ab)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                            x(b,a,c,i,j,k) = -x(a,b,c,i,j,k);
                                            x(b,a,c,k,i,j) = -x(a,b,c,i,j,k);
                                            x(b,a,c,j,k,i) = -x(a,b,c,i,j,k);
                                            x(b,a,c,i,k,j) = x(a,b,c,i,j,k);
                                            x(b,a,c,j,i,k) = x(a,b,c,i,j,k);
                                            x(b,a,c,k,j,i) = x(a,b,c,i,j,k);

                                            % (bc)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                            x(a,c,b,i,j,k) = -x(a,b,c,i,j,k);
                                            x(a,c,b,k,i,j) = -x(a,b,c,i,j,k);
                                            x(a,c,b,j,k,i) = -x(a,b,c,i,j,k);
                                            x(a,c,b,i,k,j) = x(a,b,c,i,j,k);
                                            x(a,c,b,j,i,k) = x(a,b,c,i,j,k);
                                            x(a,c,b,k,j,i) = x(a,b,c,i,j,k);

                                            % (ac)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                            x(c,b,a,i,j,k) = -x(a,b,c,i,j,k);
                                            x(c,b,a,k,i,j) = -x(a,b,c,i,j,k);
                                            x(c,b,a,j,k,i) = -x(a,b,c,i,j,k);
                                            x(c,b,a,i,k,j) = x(a,b,c,i,j,k);
                                            x(c,b,a,j,i,k) = x(a,b,c,i,j,k);
                                            x(c,b,a,k,j,i) = x(a,b,c,i,j,k);

                                            % (ac)(bc)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                            x(b,c,a,i,j,k) = x(a,b,c,i,j,k);
                                            x(b,c,a,k,i,j) = x(a,b,c,i,j,k);
                                            x(b,c,a,j,k,i) = x(a,b,c,i,j,k);
                                            x(b,c,a,i,k,j) = -x(a,b,c,i,j,k);
                                            x(b,c,a,j,i,k) = -x(a,b,c,i,j,k);
                                            x(b,c,a,k,j,i) = -x(a,b,c,i,j,k);

                                            % (ac)(ab)/[(1),(ki)(ij),(ki)(kj),(kj),(ij),(ki)]
                                            x(c,a,b,i,j,k) = x(a,b,c,i,j,k);
                                            x(c,a,b,k,i,j) = x(a,b,c,i,j,k);
                                            x(c,a,b,j,k,i) = x(a,b,c,i,j,k);
                                            x(c,a,b,i,k,j) = -x(a,b,c,i,j,k);
                                            x(c,a,b,j,i,k) = -x(a,b,c,i,j,k);
                                            x(c,a,b,k,j,i) = -x(a,b,c,i,j,k);
                                            ctQ = ctQ + 1;
                                        end

                                    end
                                end
                            end
                        end
                    end
                end
            
        end


end

function [val] = num_active(x,act_range)
    val = length(find(x >= act_range(1) & x <= act_range(end)));
end

