%% This is the processing module for DDM-UQ package
% The program will 1)transform input lag time (dt) into log space and 
%                    select 25 points. 
%                  2)process Fourier transformation for the differences 
%                    of intensity profiles on only pre-selected dt (or any
%                    user defined set), and calculate corresponding image
%                    structure function (Sqt/Dqt).
%                  3)process Fourier transformation for original intensity
%                    profile instead of the differences, which will be used 
%                    for parameter estimation in analysis module. (I_o_q_2)                  
%% Input
% I- intensity profiles, should be an array with size = L*L*len.
% option: a structure variable that may contain the following fields:
%	pxsz- the size of one pixel in unit of micron, default value is 1.
%	mindt- user determined minimum lag time, default value is 1 = mindt in
%	simulation.
%	len- number of time steps, default value is calculated by size(I,3).
%   nx- number of pixels of the image in x direction, default value is 
%       calculated by size(I,1).
%   ny- number of pixels of the image in y direction, default value is 
%       calculated by size(I,2).
%	I_o_q_2_on- parameter that control calculating I_o_q_2, which is needed
%               for estimating A(q). 1- calculate I_o_q_2, 0- doesn't
%               calculate I_o_q_2, default value is 1. 
%	sim_median_msd- median of numerical MSD value based on simulation that  
%                   comes from simulation module.
%% Output 
% Model: a structure variable that may contain the following fields:
%	dt- lag times between two consecutive images in unit of second.
%	ndt- number of frame pairs that contribute to the average/ number of
%        frame pairs that can result the same lag time.
%	selected_dt- total number of pre-selected lag times.
%	index_dt_selected- index of pre-selected lag times.
%	sz- number of pixels in square image.
%	q- wave vector in unit of um^-1.
% 	nq- number of pairs of (qx,qy) that contribute to the average/ pairs of
%       (qx,qy) that contribute to the same q where q=sqrt(qx^2+qy^2).
%	nq_index- index of pairs of (qx,qy) that contribute to the same q.
%	Sqt- also called Dqt, image structure function obtained after omitting
%        the central grid of pixels to suppress artifically large 
%        contributions due to image processing artifacts.
%	I_o_q_2- square of Fourier transformation of the original intensity
%            profile, ensamble over time t. This quantity will be used to
%            estimate A(q) in analysis module.
%   (value for input fields)

%%
function model = processing(I, varargin) 
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
if(~isfield(params,'pxsz'))
    model.pxsz = 1;
else
    model.pxsz = params.pxsz;
end

if(~isfield(params,'mindt'))
    model.mindt=1;
else
    model.mindt=params.mindt;
end 

if(~isfield(params,'len'))
    model.len = size(I,3);
else
    model.len = params.len;
end

if(~isfield(params,'nx'))
    model.nx = size(I,1);
else
    model.nx = params.nx;
end

if(~isfield(params,'ny'))
    model.ny = size(I,2);
else
    model.ny = params.ny;
end

if(~isfield(params,'I_o_q_2_on'))
    model.I_o_q_2_on = 1;
else
    model.I_o_q_2_on = params.I_o_q_2_on;
end

if(isfield(params,'sim_median_msd'))
    model.sim_median_msd = params.sim_median_msd;
end

%% Check eveything has the correct format 
if(~isnumeric(model.pxsz))
    error('the size of one pixel in unit of micron should be a numerical value');  
end

if(~isnumeric(model.mindt))
    error('the minimal lag time between 2 consecutive images should be a numerical value');  
end

if(~isnumeric(model.len))
    error('number of time steps should be a numerical value');  
end

if(~isnumeric(model.nx))
    error('number of pixels in x direction should be a numerical value');  
end

if(~isnumeric(model.ny))
    error('number of pixels in y direction should be a numerical value');  
end

if(~isnumeric(model.I_o_q_2_on))
    error('parameter that controls calculating I_o_q_2 should be a numerical value in set {0,1}');  
end

%% Generate q and dt vectors, select 25 dt points, and find corresponding sets nq and ndt
model.dt = model.mindt.*(1:model.len-1);
model.ndt = model.len-1:-1:1; 
if rem(min(model.nx, model.ny),2)==0
    model.sz = min(model.nx, model.ny)-1;
else
    model.sz = min(model.nx, model.ny);
end
edge = floor(abs(model.nx-model.ny)/2);
mid = ceil((model.sz -1)/2)+1; 
exclude = 0;

if(isfield(params,'index_dt_selected'))
    model.index_dt_selected=params.index_dt_selected;
else
% select 25 dt points and return index 
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
    index_second_part = floor(0.7*index_dt_total)+(floor(0.9*index_dt_total)-floor(0.7*index_dt_total))/5:(floor(0.9*index_dt_total)-floor(0.7*index_dt_total))/5:floor(0.9*index_dt_total);
    model.index_dt_selected = [index_dt_selected,index_second_part,index_dt_total]; 
end

model.q = (1:(model.sz-1)/2)*2*pi/(model.sz*model.pxsz); %defines q in terms of dimensional and real values
len_q = length(model.q);
[X, Y] = meshgrid(-(model.sz-1)/2:(model.sz-1)/2, -(model.sz-1)/2:(model.sz-1)/2); %make Cartesian grid
[theta ,q] = cart2pol(X, Y); %convert to polar coordinate axes
q = round(q);

model.nq = NaN(1,floor((model.sz-1)/2)); %creates a array with rows equal to half of the square dimensions
for r = 1:len_q
    model.nq(1,r) = sum(q == r,'all'); %count number of same q
end

model.nq_index = cell(1,floor((model.sz-1)/2));

for r = 1:len_q
        model.nq_index{r} = find(q == r);
end
%% Do Fourier tansformation on only pre-selected dt and calculate 1D power spectrum
Sq = double.empty;
model.Sqt = [];

count=0;
for k = model.index_dt_selected
    avg = 0;
    count = count+1;
    curr = count/length(model.index_dt_selected);
    %Tracking where the calculation is
    if rem(curr,0.2) ==0
        fprintf('Fourier transformation: %d percent \n', round(curr*100,2));
    end
    for j=1:model.ndt(k)
        imdf = I(1:model.sz,edge+1:end-edge-1, (j + k)) - I(1:model.sz,edge+1:end-edge-1, j);
        ft = fft2(imdf); %fourier transform of the difference of two images
        absft = abs(ft).^2/model.sz^2; %normalized S(qx, qy, dt) 
        avg = avg+fftshift(absft); %fftshift shifts 
    end
    avg=avg./model.ndt(k); %ensemble average
    avgRaw = avg; %avgRaw is the avg matrix before suppression
    avg(mid-exclude:mid+exclude, :) = 0;
    avg(:, mid-exclude:mid+exclude) = 0;
    avgNAN  = avg;
    avgNAN(mid-exclude:mid+exclude, :) = nan;
    avgNAN(:, mid-exclude:mid+exclude) = nan;
    
% Calculate radially averaged 1D power spectrum, S(q, dt), where q=sqrt(qx^2+qy^2)   
    for r = 1:len_q
        Sq(r,1) = nanmean( avgNAN( model.nq_index{r} ) ); 
    end
    model.Sqt = [model.Sqt, Sq];
end

%% Calculate I_o_q_2 which can be used to estimatw A(q) in analysis module
if model.I_o_q_2_on==1
    avg_I_2 = 0;
    for t = 1:model.len
        ft = fft2(I(1:model.sz,1:model.sz,t)); 
        absft = abs(ft).^2/model.sz^2; %normalized
        avg_I_2 = avg_I_2+fftshift(absft); 
    end
    avg_I_2 = avg_I_2./model.len;
    avg_I_2(mid, :) = 0;
    avg_I_2(:, mid) = 0;
    avgNAN = avg_I_2;
    avgNAN(mid, :) = nan;
    avgNAN(:, mid) = nan;
        
    model.I_o_q_2 = double.empty;
    for r = 1:len_q
       model.I_o_q_2(r,1) = nanmean( avgNAN( model.nq_index{r} ) ); 
    end
end
toc
end
