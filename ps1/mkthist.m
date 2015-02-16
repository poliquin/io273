function mkthist(fname, prices, profits, surplus, caption)
    %MKTHIST Create histograms for prices, profits, and surplus
    % Input arguments:
    %   fname    Filename for output PDF with histogram
    %   prices   Vector of prices
    %   profits  Vector of profits
    %   surplus  Vector of consumer surplus
    %   caption  Caption to display below histogram
    f = figure('PaperPosition', [.1, .2, 6, 3.5], 'PaperSize', [6.2, 4]);
    subplot(1,3,1)
    h1 = histogram(prices);
    title('Prices')
    hold on
    subplot(1,3,2)
    h2 = histogram(profits);
    title('Profits')
    subplot(1,3,3)
    h3 = histogram(surplus);
    title('Surplus')
    annotation(f, 'textbox', [0.1,0.01,0.5,0.04], ...
            'String', {caption}, ...
            'LineStyle', 'none', ...
            'FontSize', 10);
    saveas(f, fname);
end

