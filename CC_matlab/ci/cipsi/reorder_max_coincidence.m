function [det1p, det2p, occ, unocc, sign, excit_rank] = reorder_max_coincidence(det1, det2)


%     p1 = setdiff(det1,det2,'stable'); p2 = setdiff(det2,det1,'stable');
%     excit_rank = length(p1);
%     idx1 = zeros(1,excit_rank);
%     idx2 = zeros(1,excit_rank);
%     for i = 1:excit_rank
%         idx1(i) = find(det1 == p1(i));
%         idx2(i) = find(det2 == p2(i));
%     end

    [~,idx1] = setdiff(det1,det2,'stable'); occ = idx1;
    [~,idx2] = setdiff(det2,det1,'stable'); unocc = idx2;
    excit_rank = length(idx1);
    
    idx3 = setdiff(1:length(det1),idx1,'stable');
    idx4 = setdiff(1:length(det2),idx2,'stable');
    
    IDX1 = [idx1, idx3];
    IDX2 = [idx2, idx4];
    
    det1p = det1(IDX1);
    det2p = det2(IDX2);
    
    % sign is always +1???
    sign = signPermutation(IDX1)*signPermutation(IDX2);
    
        
        
    

end

