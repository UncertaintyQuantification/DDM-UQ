%% This is the analysis module for DDM-UQ package
% The program will 1) estimate noise variance, A(q,\Delta t) and B(q).
%                  2) use predictive samples from Gaussian Process 
%                     Regression (GPR) to robustly estimate the image  
%                     structure function (Dqt/Sqt), intermediate scattering
%                     function (fqt), mean squared displacement (MSD) and
%                     other quantities of interest.
%                  3) give plots for predictive Dqt and fqt for different
%                     qs with 95% predictive interval showed as shaded 
%                     area, bias detection plot for Dqt versus q to show 
%                     whether system is fully decayed or not (for not fully 
%                     decorrelated system, there's non-negligible bias for 
%                     current estimation of B(q)), unaggregated MSD for 
%                     different qs and lag times, and aggregate MSD vs lag 
%                     time with error bar that stands for 95% predictive
%                     interval. A line of true MSD for comparsion will be
%                     added if the reference value is provided. (For
%                     simulated data the generated using simulation module,
%                     median of numerical MSD will be considered as true 
%                     MSD). 
%% Input
% option: a structure variable that may contain the following fields:
%	dt- lag times between two consecutive images in unit of second
%	ndt- number of frame pairs that contribute to the average/ number of
%        frame pairs that can result the same lag time, default value is
%        length(dt).
%	q- wave vector in unit of um^-1.
%	Sqt- also called Dqt, image structure function obtained after omitting
%        the central grid of pixels to suppress artifically large 
%        contributions due to image processing artifacts.
%	subFFT- option for indicating whether the Sqt input is only a subset 
%           that is calculated using processing module, has value in {0,1}.
%           1- the given Sqt is calculated using processing module, which
%           contains data only for 25 lag times, 0- the given Sqt is the 
%           full Sqt for all possible \Delta t value. 
%	I_o_q_2- Square of Fourier transformation of the original intensity 
%            profile, ensamble over time t. 
%	index_dt_selected- index of pre-selected dt.
%   methodQ- option of q selecting method for calculating MSD, has value 
%            in {0,1}. 0- DDM-UQ paper method, also the default method; 1- 
%            the second method of selecting q, which gives a longer 
%            estimation of MSD at some circumstance. 
%   methodQ_TH- threshold paramter used in the second method of selecting
%               q, default value is 1. 
% 	sim_median_msd- median of numerical MSD value based on simulation. If
%                   this info has been provided, input data will be 
%                   considered as from a simulated experiment, and 
%                   MSD_truth = sim_median_msd.
%   MSD_truth- true MSD for used for comparsion, median of numerical MSD 
%              value from simulation or user defined reference value.
%   filename- file name for plots.
%% Output 
% Model: a structure variable that may contain the following fields:
%	A_hat- estimation of A(q,\Delta t).
%   B_hat- estimation of B(q).
%   sigma_2_0_hat- estimation of noise variance.
%   LB1- lower bound for Gaussian Process Regression parameter. 
%   LB2- lower bound for Gaussian Process Regression parameter. 
%   pred_ddm- predicted image structure function (Dqt/Sqt) for all possible 
%             time lags and first 100 qs (array with size = 100*ndt)
%   fqt- estimate intermediate scattering function (fqt) for all possible
%             time lags and first 100 qs (array with size = 100*ndt)
%   q_selected_index- index of q which are accepted for estimate MSD.
%   index_cut- contains cutting point of dt for each of the first 100 qs,
%              if have methodQ = 0.
%   selected_q_record- contains selected qs for estimation at each dt, if 
%                      have methodQ = 1.
%   t_Cut- final length of dt of MSD that we can estimate up to. 
%	median_MSD_record- median of aggregated MSD from predictive samples.
%   pre_MSD_sample- estimated MSD for different \Delta. 
%   lower_MSD_sample_95- lower bound of 95% predictive interval for MSD.
%	upper_MSD_sample_95- upper bound of 95% predictive interval for MSD.
%                              MSD.
%   (value for input fields)
%%
function model = analysis(varargin)
tic
model = struct();

if(length(varargin)>1)
    error('please put the additional inputs into a struct variable');
end

if length(varargin)==1
   params = varargin{:};
else
   params = struct();
end
%% Define default values for options
if(~isfield(params,'filename'))
    model.filename='DDM_UQ';
else
    model.filename=params.filename;
end

if(~isfield(params,'subFFT'))
    model.subFFT=1;
else
    model.subFFT=params.subFFT;
end

if(isfield(params,'dt'))
    model.dt=params.dt;
end

if(~isfield(params,'LB1'))
    model.LB1=-2E12;
else
    model.LB1=params.LB1;
end 

if(~isfield(params,'LB2'))
    model.LB2=-2E12;
else
    model.LB2=params.LB2;
end 

if(~isfield(params,'methodQ'))
    model.methodQ=0;
else
    model.methodQ=params.methodQ;
end 

if(model.methodQ==1 && ~isfield(params,'methodQ_TH'))
    model.methodQ_TH=1;
elseif(model.methodQ==1 && isfield(params,'methodQ_TH'))
    model.methodQ_TH=params.methodQ_TH;
end 

if(isfield(params,'MSD_truth'))
    model.MSD_truth=params.MSD_truth;
elseif(isfield(params,'sim_median_msd'))
    model.MSD_truth=params.sim_median_msd;
end

if(isfield(params,'index_dt_selected'))
    model.index_dt_selected=params.index_dt_selected;
else
    model.selected_dt = 25;
    log_dt = log(model.dt);
    index_dt_total = length(log_dt);
    range_log_dt_first_part = max(log_dt(1:floor(0.7*index_dt_total)))-min(log_dt);
    dist_log_dt_first_part = range_log_dt_first_part/18; 
    input_dt_benchmark = log_dt(1):dist_log_dt_first_part:log_dt(floor(0.7*index_dt_total));

    index_dt_selected = double.empty;
    i_count = 0;
    for i = 1:floor(0.7*index_dt_total)
        if(i==1|i==floor(0.7*index_dt_total)) 
            index_dt_selected = [index_dt_selected,i];
            i_count = i_count+1;
        elseif( abs(input_dt_benchmark(i_count)-log_dt(i))>=dist_log_dt_first_part )
            index_dt_selected = [index_dt_selected,i];
            i_count=i_count+1;
        end
    end
    index_second_part = floor(floor(0.7*index_dt_total)+(floor(0.9*index_dt_total)-floor(0.7*index_dt_total))/5:(floor(0.9*index_dt_total)-floor(0.7*index_dt_total))/5:floor(0.9*index_dt_total));
    model.index_dt_selected = [index_dt_selected,index_second_part,index_dt_total]; % double with dim 1*25
end

%% Check eveything has the correct format 
if(~isnumeric(model.subFFT))
    error('parameter subFFT should be a numerical value in set {0,1}');  
end

%% Input data and get estimation of noise variance, A(q,\Delta t) and B(q)
% Sqt is what people previously call Dqt
ddm_ori = params.Sqt;
q_ori = params.q;
I_o_2_ori = params.I_o_q_2;
I_o_2_ori = reshape(I_o_2_ori,[1,length(I_o_2_ori)]);

n_time = length(model.dt); 

%estimation for A, B and noise variance 
if(isfield(params,'sigma_2_0_hat'))
    model.sigma_2_0_hat = params.sigma_2_0_hat;
    model.B_hat = 2*model.sigma_2_0_hat;
else 
    if model.subFFT==0
        model.B_hat = min(mean(ddm_ori(size(ddm_ori,1),model.index_dt_selected)),min(ddm_ori(:,1))); 
    else
        model.B_hat = min(mean(ddm_ori(size(ddm_ori,1),:)),min(ddm_ori(:,1)));
    end
    model.sigma_2_0_hat = model.B_hat/2;
end

model.A_hat = 2*(I_o_2_ori - model.sigma_2_0_hat);

%% Gaussian Process Regression (GPR)
q_max_num = min(100,length(q_ori));
ddm = ddm_ori(1:q_max_num,:);
if model.subFFT==0
    sub_ddm = ddm(:,model.index_dt_selected);
else
    sub_ddm = ddm;
end
q = q_ori(1:q_max_num);
n_q = length(q);

input1 = log(q);
log_dt = log(model.dt);
input2 = log_dt(model.index_dt_selected);

p = 2;
R01 = [];
for i = 1:length(input1)
    R01(:,i)= abs(input1(i)-input1);
end
R02 = [];
for i = 1:length(input2)
    R02(:,i)= abs(input2(i)-input2);
end

n1 = length(input1);
n2 = length(input2);
num_obs = n1*n2;

output_mat = log(sub_ddm);

output = reshape(output_mat,[1,num_obs]);

param_ini = [-1,-1,0];
result_check = Neg_log_lik_eigen_with_nugget(param_ini, p, input1, input2, R01, R02, output_mat,output);

fun = @(param) Neg_log_lik_eigen_with_nugget(param, p, input1, input2, R01, R02, output_mat,output);
options = optimoptions('fmincon','GradConstr','off');  
options = optimoptions(options,'GradObj','off');
options = optimoptions(options,'MaxIterations',30);
options = optimoptions(options,'StepTolerance',10^(-5));
options = optimoptions(options,'Display','off');
options = optimoptions(options,'HessianApproximation','lbfgs');
A = [];
b = [];
Aeq = [];
beq = [];
lb = [model.LB1,model.LB2,log(0.01)];
%lb = [-inf,-inf,log(0.01)];
ub = [];
nonlcon = [];
param_opt = fmincon(fun,param_ini,A,b,Aeq,beq,lb,ub,nonlcon,options);

beta = exp(param_opt(1:p));
nu = exp(param_opt(p+1));
R1 = Matern_5_2(R01, beta(1));
R2 = Matern_5_2(R02, beta(2));
[eigen_R1_V,eigen_R1_D] = eig(R1);
[eigen_R2_V,eigen_R2_D] = eig(R2);

X = ones(num_obs,1);
q_X = size(X,2);
X_list = cell(1,q_X);
for i = 1:q_X
    X_list{i} = reshape(X(:,i),[n1,n2]);
end
U_x = [];
for i = 1:q_X   
    U_x = eigen_R1_V' * X_list{1} * eigen_R2_V;
    U_x = reshape(U_x,[1,num_obs]);
end
Lambda_tilde_inv = (1./(kron(diag(eigen_R2_D),diag(eigen_R1_D))+nu))';
Lambda_tilde_inv_U_x = Lambda_tilde_inv.*U_x;
X_R_tilde_inv_X_inv = inv(U_x * Lambda_tilde_inv_U_x');
   
output_tilde = eigen_R1_V' * output_mat * eigen_R2_V;
output_tilde = reshape(output_tilde,[1,num_obs]);
   
theta_hat = X_R_tilde_inv_X_inv * (Lambda_tilde_inv_U_x * output_tilde');

testing_input = []; 
testing_input1 = log(q);
testing_input2 = log_dt;

num_testing = n_time*n_q;  
num_testing1 = length(testing_input1);
num_testing2 = length(testing_input2);
  
testing_input(:,1) = repmat(testing_input1,1,num_testing2);
testing_input(:,2) = reshape(repmat(testing_input2,num_testing1,1),[1,num_testing]);

X_testing = ones(num_testing,1);

%Also test on a lattice 
r01 = [];
for i = 1:length(testing_input1)
    r01(:,i)= abs(testing_input1(i)-input1);
end
r02 = [];
for i = 1:length(testing_input2)
    r02(:,i)= abs(testing_input2(i)-input2);
end

r1 = Matern_5_2(r01, beta(1));
r2 = Matern_5_2(r02, beta(2));

output_mat_normalized = output-(X)'*theta_hat;
output_mat_normalized = reshape(output_mat_normalized,[n1,n2]);

output_normalize_tilde = eigen_R1_V' * output_mat_normalized * eigen_R2_V;
output_normalize_tilde = reshape(output_normalize_tilde,[1,num_obs]);
Lambda_output_normalize_tilde_mat = Lambda_tilde_inv.*output_normalize_tilde;
Lambda_output_normalize_tilde_mat = reshape(Lambda_output_normalize_tilde_mat,[n1,n2]);

v_l = r1'* eigen_R1_V;
v_u = r2'* eigen_R2_V;

predmean = v_l * Lambda_output_normalize_tilde_mat * v_u';
predmean = (X_testing)'*theta_hat+reshape(predmean,[1,num_testing]);
predmean = reshape(predmean,[num_testing1,num_testing2]); 

%Predicted Dqt and fqt 
model.pred_Dqt = exp(predmean);
model.fqt = [];
for i=1:length(q)
    model.fqt(i,:) = 1-(model.pred_Dqt(i,:)-model.B_hat)/model.A_hat(i);
end

%Variance of the mean
cstarstar = [];
for i = 1:num_testing2
    for j = 1:num_testing1
        v_here = kron(v_u(i,:),v_l(j,:));
        cstarstar((i-1)*num_testing1+j) = 1-sum(v_here.*Lambda_tilde_inv.*v_here);
    end
end 

S2 = sum(output_normalize_tilde .* Lambda_tilde_inv .* output_normalize_tilde);  
sigma2hat = S2/num_obs;

S2testing = cstarstar .* sigma2hat;
pred_var = reshape(S2testing,[num_testing1,num_testing2]);

% 95% confidence interval of the mean 
lower95_Dqt = exp(predmean+sqrt(pred_var).*norminv(0.025));
upper95_Dqt = exp(predmean+sqrt(pred_var)*norminv(0.975));

n_q_max = min(n_q,83);
n_q_min = 4;

n_sample = 300;
lower95_Dqt_sample = NaN(length(n_q_max),n_time);
upper95_Dqt_sample = NaN(length(n_q_max),n_time);
lower95_fqt_sample = NaN(length(n_q_max),n_time);
upper95_fqt_sample = NaN(length(n_q_max),n_time);

for i = 1:n_q_max 
    sample_dqt = exp(repelem(predmean(i,:),n_sample,1)+normrnd(0,1,1,n_sample)'*sqrt(pred_var(i,:)));
    sample_fqt = 1-(sample_dqt-model.B_hat)./model.A_hat(i);
    quantile_Dqt = quantile(sample_dqt,[0.025,0.975]);
    lower95_Dqt_sample(i,:) = quantile_Dqt(1,:);
    upper95_Dqt_sample(i,:)  = quantile_Dqt(2,:);
    quantile_here = quantile(sample_fqt,[0.025,0.975]);
    lower95_fqt_sample(i,:) = quantile_here(1,:);
    upper95_fqt_sample(i,:)  = quantile_here(2,:);
end

model.fqt = [];
for i=1:length(q)
    model.fqt(i,:) = 1-(model.pred_Dqt(i,:)-model.B_hat)/model.A_hat(i);
end
%% Selecting q and Sampling from predictive distribution, Nov 30 
if model.methodQ==0
    model.index_cut = zeros(1,n_q);
    index_max = zeros(1,n_q);
    % Don't accept q if it's much larger than A_hat.
    q_index_selected = [];
    l_t = round(size(model.pred_Dqt,2)*0.9);
    u_t = round(size(model.pred_Dqt,2)*0.95)-1;

    for i = n_q_min:n_q_max
        mean_end = mean(model.pred_Dqt(i,l_t:u_t));
    if mean_end<1.5*model.A_hat(i)
        q_index_selected = [q_index_selected,i];
    end
    end

    for i = 1:length(q_index_selected)
        mean_pred_Dqt_it = mean(model.pred_Dqt(q_index_selected(i),l_t:u_t));
        mean_pred_Dqt_it = min(model.A_hat(q_index_selected(i))+model.B_hat,mean_pred_Dqt_it); %should not much larger than A_hat plus B_hat
        for j = 1:n_time
            if model.pred_Dqt(q_index_selected(i),j) > 0.9*mean_pred_Dqt_it
                model.index_cut(q_index_selected(i)) = j-1;
            break
            end
        end
    end
    q_selected_index = find(model.index_cut>0);
    num_q_min = 10;
    index_t_cut = sort(model.index_cut(q_selected_index),'descend');
    
    % Sep 2023
    if length(index_t_cut)<num_q_min
         error('Accepted q is not enough, retry with "FFT_model.methodQ=1;"');
    else 
        index_t_cut = index_t_cut(num_q_min);
        index_t_truncate = find(model.index_cut(q_selected_index)>index_t_cut);
        model.index_cut(q_selected_index(index_t_truncate)) = index_t_cut;

        M_sample = 1000;
        model.pred_MSD_sample = NaN(1,max(model.index_cut));
        model.lower95_MSD_sample = NaN(1,max(model.index_cut));
        model.upper95_MSD_sample = NaN(1,max(model.index_cut));
        median_MSD_record = NaN(num_testing1,num_testing2);
        model.t_Cut = length(model.pred_MSD_sample);
        for i = 1:max(model.index_cut)
            q_index_available = find(model.index_cut>=i);
            MSD_sample_record = NaN(length(q_index_available),M_sample);
            var_record = NaN(1,length(q_index_available));
            index_here = 0;
            for j = 1:length(q_index_available)
                sample_dqt = exp(predmean(q_index_available(j),i)+sqrt(pred_var(q_index_available(j),i)).*normrnd(0,1,1,M_sample));
                index_here = index_here+1;
                index_sample_selected = find((model.A_hat(q_index_available(j))-sample_dqt+model.B_hat) > 0);
                MSD_sample_record(index_here,index_sample_selected) = (4 ./ (q(q_index_available(j)).^2)) .* log(model.A_hat(q_index_available(j))./(model.A_hat(q_index_available(j))-sample_dqt(index_sample_selected)+model.B_hat));
                var_record(index_here) = var(MSD_sample_record(index_here,index_sample_selected));
                median_MSD_record(q_index_available(j),i) = nanmedian(MSD_sample_record(index_here,index_sample_selected),'all');
            end
            MSD_sample_record_reshape = reshape(MSD_sample_record,[1,size(MSD_sample_record,1)*size(MSD_sample_record,2)]);
            quantile_here = quantile(rmmissing(MSD_sample_record_reshape),[0.025,0.975]);
            model.lower95_MSD_sample(i) = quantile_here(1);
            model.upper95_MSD_sample(i)  = quantile_here(2);
            model.pred_MSD_sample(i)  = nanmedian(MSD_sample_record,'all'); 
            if i>5
                if model.upper95_MSD_sample(i)-model.lower95_MSD_sample(i)>1.5*model.pred_MSD_sample(i)
                    model.t_Cut = i;
                    break
                end
            end
        end
    
        model.pred_MSD_sample = model.pred_MSD_sample(1:model.t_Cut);
        model.pred_MSD_sample(model.pred_MSD_sample<0) = NaN;
        model.upper95_MSD_sample = model.upper95_MSD_sample(1:model.t_Cut);
        model.upper95_MSD_sample(model.upper95_MSD_sample<0) = NaN;
        model.lower95_MSD_sample = model.lower95_MSD_sample(1:model.t_Cut);
        model.lower95_MSD_sample(model.lower95_MSD_sample<0) = NaN;
    end     
else
    t_end = floor((n_time-1)*2/3);
    model.selected_q_record = NaN(length(q_ori),t_end);
    for i = 1:t_end   
        value1 = model.A_hat(1:q_max_num)'-model.pred_Dqt(:,i)+model.B_hat;
        value2 = value1./model.pred_Dqt(:,i);
        q_index_selected = find(value2 > model.methodQ_TH & value1 >0);
        model.selected_q_record(1:length(q_index_selected),i) = q_index_selected;
    end
    t_end_cut = NaN(1,nanmax(model.selected_q_record(:,1)));
    for i = 1:nanmax(model.selected_q_record(:,1))
        ISF_i = (model.A_hat(i)-model.pred_Dqt(i,1:t_end)+model.B_hat)./model.A_hat(i);
        if isempty(min(find(ISF_i<0.2)))==0
            t_end_cut(i) = min(find(ISF_i<0.2));
        else
            t_end_cut(i) = round(0.8*t_end);
        end
        if isempty(find(islocalmin(ISF_i)))==0
            t_end_cut(i) = min(t_end_cut(i),min(find(islocalmin(ISF_i))));
        else
            t_end_cut(i) = t_end_cut(i);
        end
        model.selected_q_record(i,t_end_cut(i)+1:end) = NaN;
    end
    model.selected_q_record(1,:) = NaN;
    q23_cut = t_end_cut(4);
    for i=1:length(model.selected_q_record(1,:))
        if length(find(~isnan(model.selected_q_record(:,i))))==2
            if sum(find(~isnan(model.selected_q_record(:,i)))==[2;3])==2
                q23_cut=i;
                break
            end
        elseif length(find(~isnan(model.selected_q_record(:,i))))==1
            if find(~isnan(model.selected_q_record(:,i)))==2||find(~isnan(model.selected_q_record(:,i)))==3
                q23_cut=i;
                break
            end
        end
    end
    model.selected_q_record(2:3,1:q23_cut-1) = NaN;    
    
    M_sample = 1000;
    model.pred_MSD_sample = NaN(1,max(t_end_cut));
    model.lower95_MSD_sample = NaN(1,max(t_end_cut));
    model.upper95_MSD_sample = NaN(1,max(t_end_cut));
    median_MSD_record = NaN(num_testing1,num_testing2);
    model.t_Cut = min(max(t_end_cut(2:end)),t_end);
    for i = 1:min(max(t_end_cut(2:end)),t_end)
%         MSD_sample_record = NaN(nanmax(model.selected_q_record(:,i)),M_sample);
        MSD_sample_record = NaN(length(find(~isnan(model.selected_q_record(:,i)))),M_sample);
        for j = find(~isnan(model.selected_q_record(:,i)))'
            sample_dqt = exp(predmean(j,i)+sqrt(pred_var(j,i)).*normrnd(0,1,1,M_sample));
            index_sample_selected = find((model.A_hat(j)-sample_dqt+model.B_hat)>0);
            MSD_sample_record(j,index_sample_selected) = (4 ./ (q_ori(j).^2)) .* log(model.A_hat(j)./(model.A_hat(j)-sample_dqt(index_sample_selected)+model.B_hat));
            median_MSD_record(j,i) = nanmedian(MSD_sample_record(j,index_sample_selected),'all');
        end
        MSD_sample_record_reshape = reshape(MSD_sample_record,[1,size(MSD_sample_record,1)*size(MSD_sample_record,2)]);
        quantile_here = quantile(rmmissing(MSD_sample_record_reshape),[0.025,0.975]);
        model.lower95_MSD_sample(i) = quantile_here(1);
        model.upper95_MSD_sample(i)  = quantile_here(2);
        model.pred_MSD_sample(i)  = nanmedian(MSD_sample_record,'all');
        if i>5
            if (model.upper95_MSD_sample(i)-model.lower95_MSD_sample(i)>1.5*model.pred_MSD_sample(i)) || model.lower95_MSD_sample(i)==0
                model.t_Cut = i;
                break
            end
        end
    end
    model.pred_MSD_sample = model.pred_MSD_sample(1:model.t_Cut);
    model.pred_MSD_sample(model.pred_MSD_sample<0) = NaN;
    model.upper95_MSD_sample = model.upper95_MSD_sample(1:model.t_Cut);
    model.upper95_MSD_sample(model.upper95_MSD_sample<0) = NaN;
    model.lower95_MSD_sample = model.lower95_MSD_sample(1:model.t_Cut);
    model.lower95_MSD_sample(model.lower95_MSD_sample<0) = NaN;

end 
%% Plot 

figure ('Name',char(model.filename), 'NumberTitle','off', ...
        'Position', [400 400 1300 800])
    h=gcf;
    set(0,'DefaultAxesFontSize',12);
    set(findall(gcf,'type','text'),'FontSize',14, 'FontWeight','bold', 'fontName', 'Arial')
    set(0,'defaultaxeslinewidth','default')
    set(0,'DefaultLineLineWidth',1.5)
    markSize = 10;

subplot(2,3,1) %plotting MSD vs delta t for each q 
max_plot = round(q_max_num*0.6);
min_pred_ddm = min(min(model.pred_Dqt(round(max_plot/6):round(max_plot/6):max_plot,:)));
max_pred_ddm = max(max(model.pred_Dqt(round(max_plot/6):round(max_plot/6):max_plot,:)));
count = 1;
for i = round(max_plot/6):round(max_plot/6):max_plot
    loglog(model.dt,model.pred_Dqt(i,:),'-','Display',strcat('i = ', num2str(i)))
    hold on
    legendInfo{count} = ['q = ' num2str(round(q(i),2))]; 
    hold on
    count = count+1;
end
lgd=legend(legendInfo,'AutoUpdate','off','Location','northwest','FontSize',10);
lgd.NumColumns = 2;
for i = round(max_plot/6):round(max_plot/6):max_plot
    loglog(model.dt(model.index_dt_selected), sub_ddm(i,:), '.k', 'MarkerFaceColor','k','MarkerSize', markSize)
    hold on
end
for i = round(max_plot/6):round(max_plot/6):max_plot
    a = patch([model.dt fliplr(model.dt)], [lower95_Dqt_sample(i,:) fliplr(upper95_Dqt_sample(i,:))], [0.8 0.8 0.8],'LineStyle','none');
    uistack(a,'bottom')
    hold on
end
set(gca, 'XScale', 'log', 'YScale','log')
set(lgd,'Box','off')
xlabel('\Delta t')
ylabel('D(q,\Deltat)')
xlim([model.dt(1)-0.2*model.dt(1) model.dt(end)+0.3*model.dt(end)])
ylim([min_pred_ddm-0.2*min_pred_ddm max_pred_ddm+0.5*max_pred_ddm])
title('Predictive image structure function')
    
subplot(2,3,2)
min_pred_fqt = min(min(model.fqt(round(max_plot/6):round(max_plot/6):max_plot,:)));
max_pred_fqt = max(max(model.fqt(round(max_plot/6):round(max_plot/6):max_plot,:)));
count = 1;
for i = round(max_plot/6):round(max_plot/6):max_plot
    semilogx(model.dt,model.fqt(i,:),'-','Display',strcat('i = ', num2str(i)))
    hold on
    legendInfo{count} = ['q = ' num2str(round(q(i),2))]; 
    hold on
    count = count + 1;
end
lgd=legend(legendInfo,'AutoUpdate','off','Location','northeast','FontSize',10);
lgd.NumColumns = 2;
for i = round(max_plot/6):round(max_plot/6):max_plot
    a = patch([model.dt fliplr(model.dt)], [lower95_fqt_sample(i,:) fliplr(upper95_fqt_sample(i,:))], [0.8 0.8 0.8],'LineStyle','none');
    uistack(a,'bottom')
    hold on
end
set(gca, 'XScale', 'log')
set(lgd,'Box','off')
xlabel('\Delta t')
ylabel('f(q,\Deltat)')
xlim([model.dt(1)-0.2*model.dt(1) model.dt(end)+0.3*model.dt(end)])
ylim([min_pred_fqt-0.2*abs(max_pred_fqt) max_pred_fqt+0.3*abs(max_pred_fqt)])
title('Predictive intermediate scattering function')

subplot(2,3,3)
max_Dqt1 = max(params.Sqt(:,1));
plot(params.q, params.Sqt(:,1));
yline(0, 'r--');
xlabel('q')
ylabel('D(q,\Deltat)')
ylim([0-0.05*max_Dqt1 max_Dqt1+0.3*max_Dqt1])
title('Bias detection')

subplot(2,3,4)
min_median_MSD = min(min(median_MSD_record));
max_median_MSD = max(max(median_MSD_record));
if model.methodQ==1
    for i = 1:num_testing1
        loglog(model.dt,median_MSD_record(i,:),'-','Display',strcat('i = ', num2str(i)))
        hold on
    end
else
    for i = 1:num_testing1
        loglog(model.dt,median_MSD_record(i,:),'-','Display',strcat('i = ', num2str(i)))
        hold on
    end
end
xlabel('\Delta t')
ylabel('median MSD(q,\Deltat)')
xlim([model.dt(1)-0.2*model.dt(1) model.dt(model.t_Cut)+0.3*model.dt(model.t_Cut)])
ylim([min_median_MSD-0.2*min_median_MSD  max_median_MSD+0.5*max_median_MSD])
title('Mean squared displacement for different qs')

if(isfield(model,'MSD_truth'))
    subplot(2,3,5)
    yneg = model.pred_MSD_sample - model.lower95_MSD_sample;
    ypos = model.upper95_MSD_sample - model.pred_MSD_sample;
    errorbar(model.dt(1:model.t_Cut),model.pred_MSD_sample, yneg, ypos,'linestyle','none','Marker', '.','MarkerSize',11)
    xlabel('\Delta t')
    ylabel('MSD')
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    hold on 
    loglog(model.dt(1:model.t_Cut),model.MSD_truth(1:model.t_Cut));
    xlim([model.dt(1)-0.2*model.dt(1) model.dt(model.t_Cut)+0.2*model.dt(model.t_Cut)])
    title('Mean squared displacement')
else
    subplot(2,3,5)
    yneg = model.pred_MSD_sample - model.lower95_MSD_sample;
    ypos = model.upper95_MSD_sample - model.pred_MSD_sample;
    errorbar(model.dt(1:model.t_Cut),model.pred_MSD_sample, yneg, ypos,'linestyle','none','Marker', '.','MarkerSize',11)
    xlabel('\Delta t')
    ylabel('MSD')
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([model.dt(1)-0.2*model.dt(1) model.dt(model.t_Cut)+0.2*model.dt(model.t_Cut)])
    title('Mean squared displacement')
end
toc
end