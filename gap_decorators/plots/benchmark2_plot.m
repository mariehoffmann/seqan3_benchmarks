function foo = benchmark2_plot()
    
    M1 = csvread('benchmark2_GAP0_N8_R1024_P18.csv', 3, 0); % ok
    M2 = csvread('benchmark2_GAP1_N8_R1024_P18.csv', 3, 0); % ok
    
    subplot(1,2,1)
    title_str = 'Ins/Del at random position into UNGAPPED sequence';
    plot_one(M1, title_str, 3);
    
    title_str = 'Ins/Del at random position into GAPPED sequence';
    subplot(1,2,2)
    plot_one(M2, title_str, 3);    
    
    set(gcf,'unit','norm','position',[0 0 1 1]);
    set(gcf,'PaperPositionMode','auto');
    print('benchmark2_plot', '-dpng', '-r300');
end