function [] = test_HBar(HBar)

    addpath(genpath('/Users/harellab/Desktop/CC_matlab_tests'));
    
    load h2o-631g-HBar.mat
    %load h2o-pvdz-HBar.mat
    %load h2-pvdz-HBar.mat

    if any(any(HBar{1}{1,1} - F_mi > eps))
        fprintf('Testing Hoo... FAILED\n')
    else
        fprintf('Testing Hoo... PASSED\n')
    end

    if any(any(HBar{1}{2,2} - F_ae > eps))
        fprintf('Testing Hvv... FAILED\n')
    else
        fprintf('Testing Hvv... PASSED\n')
    end

    if any(any(HBar{1}{1,2} - F_me > eps))
        fprintf('Testing Hov... FAILED\n')
    else
        fprintf('Testing Hov... PASSED\n')
    end

    if any(any(any(any(HBar{2}{2,2,2,2} - W_abef > eps))))
        fprintf('Testing Hvvvv... FAILED\n')
    else
        fprintf('Testing Hvvvv... PASSED\n')
    end

    if any(any(any(any(HBar{2}{2,2,2,1} - W_abei > eps))))
        fprintf('Testing Hvvvo... FAILED\n')
    else
        fprintf('Testing Hvvvo... PASSED\n')
    end

    if any(any(any(any(HBar{2}{2,1,2,2} - W_amef > eps))))
        fprintf('Testing Hvovv... FAILED\n')
    else
        fprintf('Testing Hvovv... PASSED\n')
    end

    if any(any(any(any(HBar{2}{1,2,2,1} - W_mbej > eps))))
        fprintf('Testing Hovvo... FAILED\n')
    else
        fprintf('Testing Hovvo... PASSED\n')
    end

    if any(any(any(any(HBar{2}{1,2,1,1} - W_mbij > eps))))
        fprintf('Testing Hovoo... FAILED\n')
    else
        fprintf('Testing Hovoo... PASSED\n')
    end

    if any(any(any(any(HBar{2}{1,1,1,2} - W_mnie > eps))))
        fprintf('Testing Hooov... FAILED\n')
    else
        fprintf('Testing Hooov... PASSED\n')
    end

    if any(any(any(any(HBar{2}{1,1,1,1} - W_mnij > eps))))
        fprintf('Testing Hoooo... FAILED\n')
    else
        fprintf('Testing Hoooo... PASSED\n')
    end

    if any(any(any(any(any(any(HBar{3}{2,2,1,1,2,2} - W_abmief > eps))))))
        fprintf('Testing Hvvoovv... FAILED\n')
    else
        fprintf('Testing Hvvoovv... PASSED\n')
    end

    if any(any(any(any(any(any(HBar{3}{2,2,1,1,1,2} - W_abmije > eps))))))
        fprintf('Testing Hvvooov... FAILED\n')
    else
        fprintf('Testing Hvvooov... PASSED\n')
    end

    if any(any(any(any(any(any(HBar{3}{2,1,1,1,1,2} - W_amnije > eps))))))
        fprintf('Testing Hvoooov... FAILED\n')
    else
        fprintf('Testing Hvoooov... PASSED\n')
    end
    
end

