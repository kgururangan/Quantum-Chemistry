function [] = plotU(x,U,E)
    figure(2)
    plot(x,U,x,E*ones(1,length(x)))
    xlabel('x / bohr')
    ylabel('energy / Ha')
    grid on
    set(gca,'FontSize',17,'Linewidth',2,'Box','off')
    grid on
end