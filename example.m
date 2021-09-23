addpath('functions/');
%% Example for Brownian Motion
clear
option.filename = 'BM01';
option.SP = 'BM';
option.noise = 'Gaussian'; 
% simulating the process, about a few seconds
sim_model = simulation(option) 
% processing the image intensity by Fourier transformation, about 1-2 mins
FFT_model = processing(sim_model.I,sim_model) %sim_model.I-image intensity  
% compute mean squared displacement and intermediate scatter function
% about a few seconds 
GP_model = analysis(FFT_model) %based on Gaussian process regression

%% Example for experimental images
clear 
% Loading data and (Note: change path for your tiff file)
tiff_info = imfinfo('../data/30wt_30x_1.tif'); 
for k = 1:length(tiff_info)
    intensity(:,:,k) = double(imread('../data/30wt_30x_1.tif', k));
end
% giving parameter values of minimum lag time and pixel size 
option.mindt = 1/10;  %lag time
option.pxsz = 0.09697; %pixel size
FFT_model = processing(intensity,option) %it takes around 2 mins
GP_model = analysis(FFT_model)

%% Example for imported Sqt, dt and q
%% Some previous data were already fourier transformed
%% Data for an O-U process with drift
clear 
Sqt = dlmread('data/test48_Dqt.csv'); 
Sqt(1,:) = []; 
dt = dlmread('data/test48_t.csv');
q = dlmread('data/test48_q.csv');
q(1)=[];
I_o_q_2 = dlmread('data/test48_I_o_q_2.csv'); 
ndt = length(dt):-1:1;
option.Sqt = Sqt;
option.dt = dt;
option.q = q;
option.I_o_q_2 = I_o_q_2;
option.ndt = ndt;
option.subFFT = 0;
option.filename = 'test_48';
option.MSD_truth = 4*5^2*(1-0.95.^dt)+0.02^2.*dt.^2; %truth for test48
GP_model = analysis(option)
