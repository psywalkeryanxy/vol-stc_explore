clear all;close all;

%% n1001
load('tableN1001.mat');

%% exp2
load('tableexp2.mat');

yname={'pswitch';'lr';'perfCorrected';'RT(switch)'};

%pswitch,lr,perf_corrected
for whichY=1
    
    data=tableN1001;
    if whichY==1
        data.y=data.PercentSwitch;
    end
    if whichY==2
        data.y=data.LR;
    end
    if whichY==3
        data.y=data.perf_corrected;
    end
    if whichY==4
        data.y=data.rtswitch/1000;
    end
    data.x = data.log_vol_stc;
    data.x2 = (data.x).^2;
    x{1}=data.x;
    y{1}=data.y;
    
    % Fit mixed-effects model
    model_n1001=  fitlm(data,'y~ x+ x2');
    
    % Display results
    disp(model_n1001);
    exploration=data.Pswitch;
    % Get fixed effects coefficients and statistics
    fixedEffects = model_n1001.Coefficients;
    disp(fixedEffects);
    
    % Model predictions and residuals
    predictions = predict(model_n1001);
    residuals = exploration - predictions;
    
    %
    Rsquared.n1001= model_n1001.Rsquared.Adjusted;
    fprintf('R-squared: %.3f\n', Rsquared.n1001);
    
    
    
    
    
    
    %% LV
    data=tableexp2(tableexp2.condition==0,:);
    if whichY==1
        data.y=data.PercentSwitch;
    end
    if whichY==2
        data.y=data.LR;
    end
    if whichY==3
        data.y=data.perf_corrected;
    end
    if whichY==4
        data.y=data.rtswitch/1000;
    end
    data.x = log(data.vol_stc);
    data.x2 = (data.x).^2;
    x{3}=data.x;
    y{3}=data.y;
    % Fit mixed-effects model
    model_LV = fitlm(data,'y~ x+ x2');
    
    % Display results
    disp(model_LV);
    exploration=data.PercentSwitch;
    % Get fixed effects coefficients and statistics
    fixedEffects = model_LV.Coefficients;
    disp(fixedEffects);
    
    % Model predictions and residuals
    predictions = predict(model_LV);
    residuals = exploration - predictions;
    
    Rsquared.LV= model_LV.Rsquared.Adjusted;
    fprintf('R-squared: %.3f\n', Rsquared.LV);
    
    
    %% plotting with binned data
    % Setup figure
    scr_siz = get(0,'ScreenSize');
    figure('Position', scr_siz);
    
    % Define number of bins
    
    if whichY==3
        n_bins = 20;
        subplot(2,4,1);
        color=[0.6, 0.9, 0.6];
        mainfunc_plot_vsratio_bins_fitted_curve(x{3}, y{3},model_LV,n_bins,Rsquared.LV,color);
        title('LV');xlabel('log(v/s)');ylabel(yname{whichY,1});
        ylim([-0.1,0.3]);
        subplot(2,4,2);
        
        color=[0.4940, 0.1840, 0.5560];
        mainfunc_plot_vsratio_bins_fitted_curve(x{1}, y{1},model_n1001,n_bins,Rsquared.n1001,color);
        title('n1001');xlabel('log(v/s)');ylabel(yname{whichY,1});
        ylim([-0.1,0.3]);
       
    else
        
        
        
        n_bins = 20;
        subplot(2,4,1);
        color=[0.6, 0.9, 0.6];
        mainfunc_plot_vsratio_bins_fitted_curve(x{3}, y{3},model_LV,n_bins,Rsquared.LV,color);
        title('LV');xlabel('log(v/s)');ylabel(yname{whichY,1});
        subplot(2,4,2);
        color=[0.4940, 0.1840, 0.5560];
        mainfunc_plot_vsratio_bins_fitted_curve(x{1}, y{1},model_n1001,n_bins,Rsquared.n1001,color);
        title('n1001');xlabel('log(v/s)');ylabel(yname{whichY,1});
   
        
    end
    saveas(gcf,sprintf('nonlinear_%s_vsratio_%s.pdf',yname{whichY,1},date));
    close all;
    
end
