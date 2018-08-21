function foo = plot_one(M, title_str, plot_method)
    linewidth = 1.4;
    linestyles = {'--'; ':'; '-.'};
    % Vermillion Orange; SkyBlue; Green
    linecolors = [252 99 54; 0 106 224; 0 106 23]./255; 
    if (plot_method == 1)
        errorbar(log2(M(:,1)), M(:,2), M(:,3), 'LineWidth', linewidth), hold on,
        errorbar(log2(M(:,1)), M(:,4), M(:,5), linestyles{1}, 'LineWidth', linewidth), hold on,
        errorbar(log2(M(:,1)), M(:,6), M(:,7), linestyles{2}, 'LineWidth', linewidth), hold on,  
        errorbar(log2(M(:,1)), M(:,8), M(:,9), linestyles{3}, 'LineWidth', linewidth), hold on,
    elseif (plot_method == 2)
        plot(log2(M(:,1)), M(:,2), 'LineWidth', linewidth), hold on,
        plot(log2(M(:,1)), M(:,4), linestyles{1}, 'LineWidth', linewidth), hold on,
        plot(log2(M(:,1)), M(:,6), linestyles{2}, 'LineWidth', linewidth), hold on, 
        plot(log2(M(:,1)), M(:,8), linestyles{4}, 'LineWidth', linewidth), hold on, 
    else 
        semilogy(log2(M(:,1)), M(:,2), '-k', 'LineWidth', linewidth), hold on,
        semilogy(log2(M(:,1)), M(:,4), linestyles{1}, 'LineWidth', linewidth, 'Color', linecolors(1,:)), hold on,
        semilogy(log2(M(:,1)), M(:,6), linestyles{2}, 'LineWidth', linewidth, 'Color', linecolors(2,:)), hold on, 
        semilogy(log2(M(:,1)), M(:,8), linestyles{3}, 'LineWidth', linewidth, 'Color', linecolors(3,:)), hold on,         
    end
    hold off
    %set(gca,'XScale','log');
    xticks = arrayfun(@(e) pow2(e), 1:18);
    set(gca, 'XTick', xticks)
    set(gca, 'XTickLabel',[]) 
    xt = get(gca, 'XTick');
    yl = get(gca, 'YLim');
    str = cellstr( num2str(xt(:),'2^{%d}') );      %# format x-ticks as 2^{xx}
    hTxt = text(xt, yl(ones(size(xt))), str, ...   %# create text at same locations
    'Interpreter','tex', ...                   %# specify tex interpreter
    'VerticalAlignment','top', ...             %# v-align to be underneath
    'HorizontalAlignment','center');   
    title(title_str, 'FontSize', 16);
    lgd = legend('avg Gapped Sequence', 'avg Gap Vector', 'avg Anchor List', 'avg Anchor Set', ...
        'Location', 'north');
    lgd.FontSize = 16;

    ylabel('Runtime in Nanoseconds', 'FontSize', 16);
    xlabel('Sequence Length', 'FontSize', 16);
end