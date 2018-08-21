function foo = benchmark1_plot()
    
    M1 = csvread('./results/benchmark1_GAP0_N32_R1024_P18.csv', 3, 0); % ok
    M2 = csvread('./results/benchmark1_GAP1_N32_R1024_P18.csv', 3, 0); % ok
    
    ax(1) = subplot(1,2,1);
    title_str = 'Read-only on UNGAPPED sequence';
    plot_one(M1, title_str, 0);
    
    title_str = 'Read-only on GAPPED sequence';
    ax(2) = subplot(1,2,2);
    plot_one(M2, title_str, 0);    
   
    set(gcf,'unit','norm','position',[0 0 1 1]);
    set(gcf,'PaperPositionMode','auto');
    print('benchmark1_plot', '-dpng', '-r300');
end
