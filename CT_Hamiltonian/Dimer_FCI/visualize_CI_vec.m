function [] = visualize_CI_vec(dets, dets_char, pause_time)

    ylow = 1;
    yhigh = 3;
    xleft = 1;
    xright = 1.1;
    xd = 0.02;
    
    p1 = [xleft,ylow];
    p2 = [xleft+xd,ylow];
    p3 = [xright,ylow];
    p4 = [xright+xd,ylow];
    p5 = [xleft,yhigh];
    p6 = [xleft+xd,yhigh];
    p7 = [xright,yhigh];
    p8 = [xright+xd,yhigh];
    

    figure(1)
    for i = 1:size(dets,1)
        D = dets(i,:);
        for k = 1:length(D)
            hold on
            if D(k) == 1
                plot(p1(1),p1(2),'bo','MarkerSize',10,'MarkerFaceColor',[0,0,1])
                axis([0.95,1.2,ylow-0.2,yhigh+0.2])
            elseif D(k) == 2
                plot(p2(1),p2(2),'ro','MarkerSize',10,'MarkerFaceColor',[1,0,0])
                axis([0.95,1.2,ylow-0.2,yhigh+0.2])
            elseif D(k) == 3
                plot(p3(1),p3(2),'bo','MarkerSize',10,'MarkerFaceColor',[0,0,1])
                axis([0.95,1.2,ylow-0.2,yhigh+0.2])
            elseif D(k) == 4
                plot(p4(1),p4(2),'ro','MarkerSize',10,'MarkerFaceColor',[1,0,0])
                 axis([0.95,1.2,ylow-0.2,yhigh+0.2])
            elseif D(k) == 5
                plot(p5(1),p5(2),'bo','MarkerSize',10,'MarkerFaceColor',[0,0,1])
                axis([0.95,1.2,ylow-0.2,yhigh+0.2])
            elseif D(k) == 6
                plot(p6(1),p6(2),'ro','MarkerSize',10,'MarkerFaceColor',[1,0,0])
                axis([0.95,1.2,ylow-0.2,yhigh+0.2])
            elseif D(k) == 7
                plot(p7(1),p7(2),'bo','MarkerSize',10,'MarkerFaceColor',[0,0,1])
                axis([0.95,1.2,ylow-0.2,yhigh+0.2])
            else
                plot(p8(1),p8(2),'ro','MarkerSize',10,'MarkerFaceColor',[1,0,0])
                axis([0.95,1.2,ylow-0.2,yhigh+0.2])
            end
        end
        title(sprintf('Configuration %d',i))
        ll = legend(sprintf('%s',char_type(dets_char(i)))); set(ll,'Fontsize',13,'Location','NorthEast');
        pause(pause_time)
        hold off
        clf
    end

end

function [s] = char_type(n)
    if n == 0
        s = 'Other';
    elseif n == 1
        s = 'Singlet';
    elseif n == 2
        s = 'CT';
    elseif n == 3
        s = 'Triplet';
    else
        s = 'Other';
    end
end

