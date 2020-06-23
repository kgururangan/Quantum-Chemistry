function [psi_pl] = plotWF(x,psi,U,energy,ind)

    clf
    % Plotting
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.15, 0.15, 0.6, 0.6]);

    clr = jet(length(ind))*0.9; 
    %clr = lines(length(ind));
    h = zeros(1,length(ind)); leg = cell(1,length(ind));
    hold on
    for i = 1:length(ind)
        plot(x,U,'k-','Linewidth',2)
        ipl = find(U < energy(ind(i)));
        plot(x(ipl),energy(ind(i))*ones(1,length(ipl)),'color',[0.4,0,0.4],'Linewidth',2);

        psi_pl = (psi(ind(i),:)-min(psi(ind(i),:)))./(max(psi(ind(i),:))-min(psi(ind(i),:))) - 0.5;
        if length(energy) > 1
            psi_pl = psi_pl*mean(diff(energy))*0.7+abs(energy(ind(i)));
        else
            psi_pl = psi_pl*0.2*U(end)+abs(energy(ind(i)));
        end
     
        h(i) = plot(x,psi_pl,'color',clr(i,:),'Linewidth',2);
        leg{i} = sprintf('E_{%d} = %4.2f cm^{-1}',ind(i),energy(ind(i)));
        xlabel('$x - x_{eq}$ / $\AA$','Interpreter','latex')
        ylabel('$U(x)$ / cm$^{-1}$','Interpreter','latex')
        set(gca,'FontSize',17,'Linewidth',2,'Box','off');
        axis tight
    end
    hold off
    ll = legend(h,leg); set(ll,'FontSize',10,'Location','eastoutside');
    
end

