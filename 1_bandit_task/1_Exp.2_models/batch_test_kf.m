clear all;close all;


load('alldatabandit.mat');
%LV then HV
%[LV;HV];

load('lap_kf_v6p25.mat');

allsub_sigma=[];
allsub_omega = [];
allsubunc = [];
allsublr = [];
allsubval_1=[];
allsubval_2=[];
allsubval_3=[];

allsub_pchoice_1 = [];
allsub_pchoice_2 = [];
allsub_pchoice_3 = [];

for whichsub=1:length(alldata)
    
    
    data = alldata{whichsub,1};
    
    ux       = @(x)(1/(1+exp(-x)));
    
    params = cbm.output.parameters(whichsub,:);
    
    
    lux      = @(x)(1/(1+exp(-x)));
    sigma    = lux(params(1));
    omega    = lux(params(2));
    params_response = params(3);
    
    choice   = data.choice;
    outcome  = data.outcome;
    
    
    
    [dv,unc,lr]  = kf(outcome,sigma,omega);
    [pchoice] = fit_p(dv,params_response );
    %dv=2*dv-1;
    allsubunc = [ allsubunc;unc(:,1)'];
    allsubval_1=[allsubval_1;dv(:,1)'];
    allsubval_2=[allsubval_2;dv(:,2)'];
    allsubval_3=[allsubval_3;dv(:,3)'];
    
    allsublr = [allsublr;lr(:,1)'];
    
    allsub_sigma=[allsub_sigma;sigma];
    allsub_omega = [allsub_omega;omega];
    
    
    allsub_pchoice_1 = [allsub_pchoice_1;pchoice(:,1)'];
    allsub_pchoice_2 = [allsub_pchoice_2;pchoice(:,2)'];
    allsub_pchoice_3 = [allsub_pchoice_3;pchoice(:,3)'];
end

save results_KF


%save results_kf_v6p25 all*

function [dv,walltrial,kalltrial] = kf(y,sigma,omega)

[nt,nq] = size(y);

m       = zeros(1,nq)+0.33;
w       = sigma*ones(1,nq);
dv      = nan(nt,nq);
walltrial=[];kalltrial=[];
for t  = 1:nt
    dv(t,:)     = (m(1,:));
    
    k           = (w+sigma)./(w+sigma + omega);
    delta       = y(t,:) - m;
    m           = m + k.*delta;
    w           = (1-k).*(w+sigma);
    walltrial=[walltrial;w];
    kalltrial=[kalltrial;k];
end

end

% Y        = choice==1;
% loglik   = fit_A_response(dv,Y,params_response);
function [p] = fit_p(X,params)
beta     = exp(params(1));

sumX=sum(X,2);
adaptX = X./sumX;

z        = adaptX*beta;
f        = (1./(1+exp(-z)));

p=f;

end

