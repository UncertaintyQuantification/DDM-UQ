%% This is the simulation module of DDM-UQ package
% The program will 1)simulate 2-D particle movement. 
%                  2)calculate the numerical MSD.
%                  3)save particle trajectories, parameter values and 
%                    intensity profiles in current directory (if save=1).
%                  4)give plots for simulated particle trajectories and
%                    numerical MSD.
%% Input
% option: a structure variable that may contain the following fields:
%   filename- user defined file name for current simulation, default value
%             is 'sim01'.
%	save- option for saving, has value in {0,1}. 1- save all parameter 
%         values in a text file; particle trajectories, and intensity 
%         profiles as .mat file in current directory. 0- results will not 
%         be saved in seperate files. Default value is 0.
%	N- number of particles, default value is 50. 
%	L- side length of the image, assuming square image with L*L pixels,  
%      default value is 480 (480*480 pixels).
%	len- number of time steps, default value is 1000.
%	SP- type of stochastic process used for simulation, default is 'BM'
%	    (Brownian motion/simple diffusion), another option is 'OU'
%       (Ornstein–Uhlenbeck process).
%	step- distance moved per timestep = sqrt(2D*timestep), default value
%         is 2.
%	drift- option for adding a drift, has value in {0,1}, 1- with a drift, 
%          0- without a drift, default value is 0.       
%	drift_dir- option for direction of the drift, has value in {0,1}. 
%              1- direction of the drifts are distinct to each 
%              particle, fixed in direction theta_i for particle i. 
%              0- same drift for all particles. Default value is 0.
%	mu- mean velocity of the drift, default value is 0.1.
%	rho- correlation between successive step and previous step in O-U 
%        process, default value is 0.95. 
%	samepos- option for generating all particles in same location, has 
%            value in {0,1}. 1- initial positions for all N particles are 
%            the same (optically densen sample). 0- initial positions are 
%            random (optically dilute sample). Default value is 0.
%	pos0- initial [x,y] positions for all N particles. If samepos=1,  
%         select a random location that is in the middle of the frame
%         [4/10-6/10]L*[4/10-6/10]L, and use that for all particles. If 
%         samepos=0, randomly generate N locations within L*L frame. Or can
%         be any user defined value (need to be a N*2 array).
%	noise- option for adding time varying background noise, has value in 
%          {'uniform','Gaussian'}, default value is 'uniform'.
%	sigma- width of the Gaussian, default value is 2.
%	Imax- maximum intensity at the center of the particle (at different 
%         time steps), default value is 255. User defined input for this 
%         parameter can be a single number or a vector with length equals
%         to len if time varying signal is perfered. 
%	I0- background intensity (at different time steps), default value is 
%       20. User defined input for this parameter can be a single number or
%       a vector with length equals to len if time varying noise is 
%       perfered.
%% Output 
% Model: a structure variable that may contain the following fields:
%	pos- particle positions for all N particles in all time steps, an array 
%        with size = N*len*2.
%	I- intensity profiles, an array with size = L*L*len.
%   sim_median_msd- median of numerical MSD value based on simulation. This
%                   info can be used in analysis module for comparsion. 
%   (value for input fields)
% If save = 1:
%   1 text file that contains all used parameter values.
%   2 .mat files- 1) particle positions for all time steps.
%                 2) intensity profiles.
%% Examples (Copy the code to example.m file to run following examples.)
% For function details and reference information, see: 
% Gu, M., Luo, Y., He, Y., Helgeson, M. E., & Valentine, M. T. (2021). 
% Uncertainty quantification and estimation in differential dynamic microscopy. 
% Physical Review E, just accepted, arXiv preprint arXiv:2105.01200.

% %Example 1: simple diffusion 
% clear
% option.filename = 'BM01'; 
% option.SP = 'BM';
% sim_model = simulation(option)

% %Example 2: diffusion with a drift (optically dilute sample)
% clear
% option.SP = 'BM';
% option.drift = 1;
% option.step = .56;
% option.mu = 0.1;
% sim_model = simulation(option)

% %Example 3: diffusion with a drift (optically dense sample)
% clear
% option.SP = 'BM';
% option.drift = 1;
% option.step = .56;
% option.mu = 0.1;
% option.samepos = 1;
% sim_model = simulation(option)

% %Example 4: O-U process with a drift 
% clear
% option.SP = 'OU';
% option.drift = 1;
% option.step = 5;
% option.mu = 0.02;
% option.rho = 0.95;
% sim_model = simulation(option)

% %Example 5: O-U process with drifts distinct to each particle  
% clear
% option.SP = 'OU';
% option.noise = 'Gaussian'; 
% option.drift = 1;
% option.drift_dir = 1;
% option.step = 5;
% option.mu = 0.02;
% option.rho = 0.95;
% option.N = 50;
% option.L = 480;
% option.len = 1000;
% option.pos0 = option.L/8+3*option.L/4.*rand(option.N,2);
% sim_model = simulation(option)
%% 
function model= simulation(varargin) 
tic
if(length(varargin)>1)
    error('please put the additional inputs into a struct variable');
end

if length(varargin)==1
   params=varargin{:};
else
    params=struct();
end

%% Define default values for options
if(~isfield(params,'filename'))
    model.filename = 'sim01';
else
    model.filename = params.filename;
end

if(~isfield(params,'save'))
    model.save = 0;
else
    model.save = params.save;
end

if(~isfield(params,'N'))
    model.N = 50;
else
    model.N = params.N;
end

if(~isfield(params,'L'))
    model.L = 480;
else
    model.L = params.L;
end

if(~isfield(params,'len'))
    model.len = 1000;
else
    model.len = params.len;
end

if(~isfield(params,'SP'))
    model.SP = 'BM';
else
    model.SP = params.SP;
end

if(~isfield(params,'step'))
    model.step = 2;
else
    model.step = params.step;
end

if(~isfield(params,'drift'))
    model.drift = 0;
else
    model.drift = params.drift;
end

if(model.drift==1)
    if(~isfield(params,'drift_dir'))
        model.drift_dir = 0;
    else
        model.drift_dir = params.drift_dir;
    end
end

if(model.drift==1)
    if(~isfield(params,'mu'))
        model.mu = 0.1;
    else
        model.mu = params.mu;
    end
elseif(model.drift==0)
    model.mu = 0;
end

if(strcmp(model.SP,'OU'))
    if(~isfield(params,'rho'))
        model.rho = 0.95;
    else
        model.rho = params.rho;
    end
end

if(~isfield(params,'samepos'))
    model.samepos=0;
else
    model.samepos = params.samepos;
end

if(~isfield(params,'pos0') && model.samepos==0)
    model.pos0 = model.L.*rand(model.N,2);
elseif(~isfield(params,'pos0') && model.samepos==1)
    model.pos0 = model.L/2+model.L/10.*repmat(rand(1,2),model.N,1);
elseif(isfield(params,'pos0'))
    model.pos0 = params.pos0;
end

if(~isfield(params,'noise')) %Gaussian or uniform
        model.noise = 'uniform';
else
        model.noise = params.noise; 
end

if(~isfield(params,'sigma'))
    model.sigma = 2;
else
    model.sigma = params.sigma;
end

if(~isfield(params,'Imax'))
    model.Imax = 255;
else
    model.Imax = params.Imax;
end

if(~isfield(params,'I0'))
    model.I0 = 20;
else
    model.I0 = params.I0;
end

%% Check eveything has the correct format
if(~ischar(model.filename))
    error('file name should be a character value');  
end

if(~isnumeric(model.save))
    error('saving option should be a numerical value in set {0,1}');  
end

if(~isnumeric(model.N))
    error('number of particles should be a numerical value');  
end

if(~isnumeric(model.L))
    error('number of pixels in one side should be a numerical value');  
end

if(~isnumeric(model.len))
    error('number of time steps should be a numerical value');  
end

if(~ischar(model.SP))
    error('type of particle movement we simulate from should be a character value in set {BM, OU}');  
end

if(~isnumeric(model.step))
    error('distance moved per timestep should be a numerical value');  
end

if(~isnumeric(model.drift))
    error('drift option should be a numerical value in set {0,1}');  
end

if(~isnumeric(model.samepos))
    error('option for whether initial positions are the same or not should be a numerical value in set {0,1}');  
end

if(~isnumeric(model.sigma))
    error('sigma should be a numerical value');  
end

if(~isnumeric(model.Imax))
    error('maximum intensity should be numerical');  
end


if(~isnumeric(model.I0))
    error('minimum intensity for the background should be numerical');  
end

%% Simulate 2D particle movement and calculate numerical MSD
% return positions for all particles at each time step
if(strcmp(model.SP,'BM'))
    model.pos = diffusion_with_drift(model.len,model.step,model.pos0,model.N,model.mu); %pos is an array with dim=N*len*2
elseif(strcmp(model.SP,'OU'))
    model.pos =  OU_diff_drift_every_2_time_points(model.len,model.N,model.pos0,model.drift_dir,model.step,model.rho,model.mu);
end

% numerical MSD
msdi = nan(model.N,model.len-1); %MSD for each particle 
for dt = 1:model.len-1
        nt = model.len - dt;
        xdiff = model.pos(:,1:nt,1)-model.pos(:,(1+dt):(nt+dt),1);
        ydiff = model.pos(:,1:nt,2)-model.pos(:,(1+dt):(nt+dt),2);
        msdi(:,dt) = nanmean(xdiff.^2 + ydiff.^2,2);
end
msd = nanmean(msdi,1);
model.sim_median_msd = nanmedian(msdi,1);

%% Construct intensity profiles 
if length(model.I0)>1
    if(strcmp(model.noise,'uniform'))
        Rand = rand(model.L,model.L,model.len)-0.5;
        model.I = reshape(model.I0.*reshape(Rand,model.L*model.L,model.len),model.L,model.L,model.len);       
    elseif(strcmp(model.noise,'Gaussian'))
        Rand = randn(model.L,model.L,model.len);
        model.I = reshape(sqrt(model.I0).*reshape(Rand,model.L*model.L,model.len),model.L,model.L,model.len);       
    end
else
    if(strcmp(model.noise,'uniform'))
        model.I = model.I0*(rand(model.L,model.L,model.len)-0.5);       
    elseif(strcmp(model.noise,'Gaussian'))
        model.I = sqrt(model.I0)*(randn(model.L,model.L,model.len));       
    end
end

Ic = model.Imax.*ones(model.N,1); %unifrom brightness at the center of the particle
I_matrix=reshape(model.I,[model.L*model.L,model.len]);
if length(model.Imax)>1
    for k=1:model.len
        for j=1:model.N       
            [index_fill, Ip_fill]=particle_intensity(model.pos(j,k,1),model.pos(j,k,2),model.sigma,Ic(j,k),model.L);    
            legitimate_index=((index_fill>0)&(index_fill<(model.L*model.L)) );
            if(size(legitimate_index,1)>0)
                I_matrix(index_fill(legitimate_index),k)=I_matrix(index_fill(legitimate_index),k)+Ip_fill(legitimate_index);
            end
        end
    end
else
    for k=1:model.len
        for j=1:model.N       
            [index_fill, Ip_fill]=particle_intensity(model.pos(j,k,1),model.pos(j,k,2),model.sigma,Ic(j),model.L);    
            legitimate_index=((index_fill>0)&(index_fill<(model.L*model.L)) );
            if(size(legitimate_index,1)>0)
                I_matrix(index_fill(legitimate_index),k)=I_matrix(index_fill(legitimate_index),k)+Ip_fill(legitimate_index);
            end
        end
    end    
end 
model.I=reshape(I_matrix,[model.L,model.L,model.len]);

%% Plot 

figure ('Name',model.filename, 'NumberTitle','off', ...
        'Position', [500 500 700 310])
    h=gcf;
    set(0,'DefaultAxesFontSize',12);
    set(findall(gcf,'type','text'),'FontSize',14, 'FontWeight','bold', 'fontName', 'Arial')
    set(0,'defaultaxeslinewidth','default')
    set(0,'DefaultLineLineWidth',1.5)
    markSize =5;

subplot(1,2,1)
for i = 1:model.N
    %figure(1)
    plot(model.pos(i,:,1),model.pos(i,:,2))
    hold on 
    pbaspect([1 1 1])
    axis ([0 model.L 0 model.L])
end
title('Particle trajectories')

subplot(1,2,2)
for i = 1:model.N
    loglog(1:model.len-1,msdi(i,:),'-','Display',strcat('i = ', num2str(i)))
    hold on
end
loglog(1:model.len-1,msd,'k-','LineWidth',2)
xlabel('\Delta t [timestep]')
ylabel('MSD [pixel^2]')
title('Mean Square Displacement')

%% Book keeping
if model.save==1
    str = pwd;
    fileID = fopen(strcat(str,'/',model.filename,'.txt'), 'w');
    fprintf(fileID,'image size = %d \n',model.L);
    fprintf(fileID,'particle number = %d \n',model.N);
    fprintf(fileID,'number of time steps taken = %d\n',model.len);
    fprintf(fileID,'pixels moved per timestep (= sqrt(2D*timestep)) = %d \n',model.step);
    fprintf(fileID,'width of the Gaussian = %d \n',model.sigma);
    fprintf(fileID,'maximum intensity = %0.2f \n', model.Imax);
    fprintf(fileID,'image noise level = %0.2f \n', model.I0);
    fclose(fileID);
    pos = model.pos;
    I = model.I;
    %save(strcat(str,'/',model.filename,'_pos.mat'),'pos','-v6') 
    %save(strcat(str,'/',model.filename,'_I.mat'),'I','-v6')
end

toc
end