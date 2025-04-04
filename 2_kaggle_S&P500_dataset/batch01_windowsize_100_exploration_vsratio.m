clear all;close all;
scr_siz = get(0,'ScreenSize');

winsize=[100:10:150]';num_windows=10;
orange_map = [linspace(1,0.8,num_windows)', linspace(0.5,0.2,num_windows)', linspace(0,0,num_windows)'];
figure('Position', [50, 50, 700, 400]);
for whichWinDowSize=1:length(winsize)
    % load data
    data = readtable('all_stocks_5yr.csv');  %
    
    % preprocessing...
    data.date = datetime(data.date);
    
    
    
    % top 50 stocks
    N = 50;
    volume_mean = grpstats(data, 'Name', 'mean', 'DataVars', 'volume');
    [~, idx] = sort(volume_mean.mean_volume, 'descend');
    top_stocks = volume_mean.Name(idx(1:N));
    
    % this top stocks'data
    data = data(ismember(data.Name, top_stocks), :);
    
    % get price matrices (time x stocks)
    %
    dates = unique(data.date);
    symbols = unique(data.Name);
    
    % close prices
    prices = nan(length(dates), length(symbols));
    for i = 1:length(symbols)
        stock_data = data(strcmp(data.Name, symbols{i}), :);
        [~, date_idx, ~] = intersect(dates, stock_data.date);
        prices(date_idx, i) = stock_data.close;  % use the close price
    end
    
    % check print
    fprintf('Loaded price data:\n');
    fprintf('Time period: %s to %s\n', ...
        datestr(dates(1)), datestr(dates(end)));
    fprintf('Selected stocks: \n');
    for i = 1:length(symbols)
        fprintf('%s\n', symbols{i});
    end
    fprintf('Number of trading days: %d\n', size(prices, 1));
    
    % data quality
    missing_data = sum(isnan(prices(:)));
    if missing_data > 0
        warning('Found %d missing values', missing_data);
        % nan
        for i = 1:size(prices, 2)
            prices(:,i) = fillmissing(prices(:,i), 'linear');
        end
    end
    
    % analyses
    
    % Parameters
    window_size = winsize(whichWinDowSize,1);
    
    % returns
    returns = (prices(2:end,:) - prices(1:end-1,:)) ./ prices(1:end-1,:);
    
    % volatility,stochasticity
    [volatility, stochasticity] = calculate_uncertainty(returns, window_size);
    
    %stochasticity=1*10-16e
    filterstc=find(stochasticity<2.2204e-15);
    volatility(filterstc)=[];
    stochasticity(filterstc)=[];
    
    ratio = log(volatility ./ stochasticity);
    
    % exploration (portfolio diversity)
    exploration = calculate_portfolio_diversity(returns, window_size);
    exploration(filterstc)=[];
    %remove outlier (1%,99% quantile)
    valid_idx = ~isnan(ratio) & ~isnan(exploration);
    ratio_clean = ratio(valid_idx);
    exploration_clean = exploration(valid_idx);
    
    %remove outlier (1%,99% quantile)
    [ratio_filtered, exp_filtered] = remove_outliers(ratio_clean, exploration_clean);
    %% delete the extreme values
    
    IDXF=find(ratio_filtered>100);
    ratio_filtered(IDXF)=[];
    exp_filtered(IDXF)=[];
    %
    datatable=table;
    datatable.y= exp_filtered;
    datatable.x=ratio_filtered;
    datatable.x2=(datatable.x).^2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lmresults=fitlm(datatable,'y~x+x2');
    
    
    
    
    %% plotting with binned data
    % Setup figure
    
    x{1}=datatable.x;
    y{1}=datatable.y;
    
    color=orange_map(whichWinDowSize,:);
    n_bins=20;
    % Plot the fitted curve with the same orange spectrum color
    
    subplot(1,3,1)
    mainfunc_plot_vsratio_bins_fitted_curve(x{1}, y{1},lmresults,n_bins,lmresults.Rsquared.Adjusted,color);
    title('S&P500 5yrs');xlabel('log(volstc ratio)');ylabel('portfolio diversity(exploration)');

    hold on
    
    
end
% Add a colorbar to show window size scale
colormap(orange_map);
c = colorbar;
c.Label.String = 'Window Size (days)';
c.Ticks = linspace(0, 1, 6);
c.TickLabels = {num2str(winsize(1)), '', '', '', '', num2str(winsize(end))};

saveas(gcf,'SP_500_winsize_100_150_explore_vsratio.pdf');


function [volatility, stochasticity] = calculate_uncertainty(returns, window_size)
% Calculate rolling volatility
volatility = movstd(returns, window_size);
volatility = mean(volatility, 2); % avg

% Calculate stochasticity using a more robust method
stochasticity = zeros(size(returns, 1), 1);
for t = window_size:length(returns)
    window_data = returns(t-window_size+1:t, :);
    % robust stoc
    stochasticity(t) = calculate_robust_stochasticity(window_data);
end

% make sure no 0
stochasticity = max(stochasticity, eps);
% make sure no extreme value
stochasticity = winsorize(stochasticity, 0.01); % 1% extremevalue winsorization
end

function stoch = calculate_robust_stochasticity(data)
% 1. get return dispersion
dispersion = std(data(:));

% 2. get autocorrelation
[acf, lags] = autocorr(mean(data,2), 'NumLags', min(5,size(data,1)-1));
acf = acf(2:end); % removelag 0

% 3. get runs test
runs_stat = calculate_runs_stat(data);

% 4. summarize
stoch = mean([dispersion, 1-mean(abs(acf)), runs_stat]);

% 
stoch = max(stoch, eps);
end

function stat = calculate_runs_stat(data)
% runs test
returns_mean = mean(data(:));
binary_seq = data(:) > returns_mean;
runs = diff([0; binary_seq; 0] ~= 0);
num_runs = sum(runs);

% standarliz
n1 = sum(binary_seq);
n2 = length(binary_seq) - n1;
exp_runs = 2*n1*n2/length(binary_seq) + 1;
stat = num_runs/exp_runs;

% 0-1
stat = max(0, min(1, stat));
end

function data_out = winsorize(data_in, threshold)
% Winsorization
lower = prctile(data_in(data_in > 0), threshold*100);
upper = prctile(data_in(data_in > 0), (1-threshold)*100);
data_out = data_in;
data_out(data_out < lower) = lower;
data_out(data_out > upper) = upper;
end



function entropy = calculate_entropy(data)
% Approximate entropy calculation
data_mean = mean(data(:));
data_std = std(data(:));
% standarize
norm_data = (data - data_mean) / data_std;
% get entropy
[counts, ~] = histcounts(norm_data, 'BinMethod', 'fd');
probs = counts / sum(counts);
probs = probs(probs > 0);
entropy = -sum(probs .* log2(probs));
end

function diversity = calculate_portfolio_diversity(returns, window_size)
% Calculate portfolio diversity (exploration measure)
[n_days, n_stocks] = size(returns);
diversity = zeros(n_days, 1);

for t = window_size:n_days
    % get weight
    window_returns = returns(t-window_size+1:t, :);
    perf = mean(window_returns, 1);
    % non-negative weight
    weights = perf - min(perf) + eps;
    weights = weights / sum(weights);
    % diversity (1 - HHI)
    diversity(t) = 1 - sum(weights.^2);
end
end

function [ratio_filtered, exp_filtered] = remove_outliers(ratio, exploration)
% Remove outliers based on percentiles
ratio_prctile = prctile(ratio, [1 99]);
valid_idx = ratio >= ratio_prctile(1) & ratio <= ratio_prctile(2);
ratio_filtered = ratio(valid_idx);
exp_filtered = exploration(valid_idx);
end

