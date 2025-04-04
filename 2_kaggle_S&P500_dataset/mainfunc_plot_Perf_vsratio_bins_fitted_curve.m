function mainfunc_plot_Perf_vsratio_bins_fitted_curve(x,y,model,n_bins, windowsize,color)
% n1001 data (purple)
[bin_centers_n1001, bin_means_n1001, bin_sems_n1001] = bin_data(x, y, n_bins);
plot(bin_centers_n1001, bin_means_n1001, 'o', 'Color', color, ...
    'MarkerSize', 1);

hold on
% n1001 fit
coef = model.Coefficients.Estimate;
pValues = model.Coefficients.pValue;

% Create smooth x range for the curve
x_smooth = linspace(min(x), max(x), length(x));
% Calculate predicted values using quadratic equation
y_pred = coef(1) + coef(2)*x_smooth + coef(3)*(x_smooth.^2);
plot(x_smooth, y_pred, '-', 'Color', color, 'LineWidth', 1);
% xlim([-3.5,-2]);
%ylim([-0.05,0.05]);

% Remove box and add minimal axes
box off;
ax = gca;
ax.XAxisLocation = 'bottom';
ax.YAxisLocation = 'left';
% Add p-values and R² to plot
% Convert p-values to stars for significance
function sig_str = get_stars(p)
    if p < 0.001
        sig_str = '***';
    elseif p < 0.01
        sig_str = '**';
    elseif p < 0.05
        sig_str = '*';
    else
        sig_str = 'n.s.';
    end
end

% Format p-values with scientific notation if very small
function p_str = format_pvalue(p)
    if p < 0.001
        p_str = sprintf('p = %.2e', p);
    else
        p_str = sprintf('p = %.3f', p);
    end
end

% Add text annotations
text(-3, -0.0002, sprintf('Linear term: %s (%s)', get_stars(pValues(2)), format_pvalue(pValues(2))), ...
    'FontSize', 8, 'FontName', 'Arial');
text(-3, -0.0004, sprintf('Quadratic term: %s (%s)', get_stars(pValues(3)), format_pvalue(pValues(3))), ...
    'FontSize', 8, 'FontName', 'Arial');
text(-3, -0.00060, sprintf('window size=%d days', windowsize), ...
    'FontSize', 8, 'FontName', 'Arial');
% Function to bin data and calculate means
    function [bin_centers, bin_means, bin_sems] = bin_data(x, y, n_bins)
        [~, edges] = histcounts(x, n_bins);
        bin_centers = movmean(edges, 2, 'Endpoints', 'discard');
        bin_means = zeros(1, n_bins-1);
        bin_sems = zeros(1, n_bins-1);
        
        for i = 1:(length(edges)-1)
            bin_idx = x >= edges(i) & x < edges(i+1);
            if i == length(edges)-1
                bin_idx = x >= edges(i) & x <= edges(i+1);
            end
            bin_means(i) = mean(y(bin_idx), 'omitnan');
            bin_sems(i) = std(y(bin_idx), 'omitnan') / sqrt(sum(bin_idx));
        end
    end
end