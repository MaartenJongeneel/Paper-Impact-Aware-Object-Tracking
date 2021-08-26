close all; clearvars; clc; 
addpath(genpath('static')); addpath('Functions'); 
set(groot,'defaulttextinterpreter','latex'); set(groot,'defaultAxesTickLabelInterpreter','latex');  set(groot,'defaultLegendInterpreter','latex');
%% Cuboid Tracking using Particle Filter with Constant Velocity model.
config = 1; %Select the configuration

if config == 1; fname = "config01"; elseif config == 2; fname = "config02"; end
load(append(fname,'/K.mat'));              %Load camera parameters
load(append(fname,'/box.mat'));            %Load the model of the box
load(append(fname,'/GT.mat'));             %Load the ground truth data
FileList = dir(fullfile(append('static/',fname,'/Test_data'), '*.png')); %Load the data set

%% Settings
Npart    = 500;      %Number of particles used   [-]
maxt     = 65;       %Run to this frame          [-]
alpha    = 0.9;      %UKF : point scaling parameter
beta     = 1;        %UKF : scaling parameter for higher order terms of Taylor series expansion
kappa    = 0.5;      %UKF : sigma point selection scaling parameter 
DoSave   = true;    %Decide if you want to save 

if config ==1
    %Constant parameters
    const.Noise = 0.001*diag([1.5 1.5 0.5 1 1 1.5]);
    impacts = [16 31 39 49 58]; %Frame at which an impact occurs
    
    %Process noise covariance and measurement noise covariance
    Qv = 5e-4*diag([15 15 50 1 10 1 10 50 1 30 30 100]); %[mm]  Process noise covariance
    Rv = 5e-4*diag([5 5 5 1 1 1]);  %[mm]  Measurement noise covariance

    %Initial state covariance 
    PR = 5e-4*diag([1 1 1]);    Po = 1e-4*diag([1 1 1]);
    Pv = 1e-7*diag([1 1 1]);    Pw = 1e-7*diag([1 1 1]);

    %Initial pose and velocity
    AR_B = GT{1}(1:3,1:3);   %Initial Rotation
    Ao_B = GT{1}(1:3,4);     %Initial position [m]
    Av_B = [ -2;    0;   0]; %Initial linear velocity [m/s]
    Aw_B = [  5;    2;   0]; %Initial angular velocity [rad/s] 
    
elseif config == 2
    %Constant parameters
    const.Noise = 0.001*diag([1.7 11.9 1 84 17 0.4]);
    impacts = [14 20 28 35 45]; %Frame at which an impact occurs
    
    %Process noise covariance and measurement noise covariance
    Qv = 4.5e-4*diag([15 15 50 1 10 5 10 50 1 50 50 100]); %[mm]  Process noise covariance
    Rv = 5e-4*diag([5 5 5 1 1 1]);  %[mm]  Measurement noise covariance

    %Initial state covariance 
    PR = 5e-4*diag([1 1 1]);    Po = 1e-4*diag([1 1 1]);
    Pv = 1e-7*diag([1 1 1]);    Pw = 1e-7*diag([1 1 1]);

    %Initial pose and velocity
    AR_B = GT{1}(1:3,1:3);   %Initial Rotation
    Ao_B = GT{1}(1:3,4);     %Initial position         [m]
    Av_B = [  1; -1.2;   0]; %Initial linear velocity  [m/s]
    Aw_B = [  5;    2;   0]; %Initial angular velocity [rad/s]
end

%% Reference histograms
HA = HIS(imread('Red.png'),12,12,4);
HB = HIS(imread('Green.png'),12,12,4);
HC = HIS(imread('Cyan.png'),12,12,4);
HD = HIS(imread('Magenta.png'),12,12,4);
HE = HIS(imread('Yellow.png'),12,12,4);
HF = HIS(imread('Blue.png'),12,12,4);
Href = {HA,HB,HC,HD,HE,HF};

%% 1: INITIALIZATION
time = tic();

%Sample inital set of particles
img=imread(fullfile(append('static/',fname,'/Test_data'), FileList(1).name)); %Load the image
hsi = rgb2hsi(img,12,12,4);                          %Convert image to HSI bins (12x12x4)
for i = 1:Npart
    xiR = zeros(3,1) + sqrtm(PR)*randn(3,1);
    xio = zeros(3,1) + sqrtm(Po)*randn(3,1);
    xiv = zeros(3,1) + sqrtm(Pv)*randn(3,1);
    xiw = zeros(3,1) + sqrtm(Pw)*randn(3,1);
    
    X{1}{1,i}  = AR_B*expm(hat(xiR));
    X{1}{2,i}  = Ao_B+xio;
    X{1}{3,i}  = Av_B+xiv;
    X{1}{4,i}  = Aw_B+xiw;
    
    P{1}{1,i}  = blkdiag(PR,Po,Pv,Pw);              
    [wK(i),~,~,~]=likelihood(X{1}(:,i),Href,hsi,K,box); 
end 

%Compute the weighted average of the inital frame
wK=wK./sum(wK);       %Normalize the weights
[~,indx] = max(wK);   %Index of the particle with highest weight
wmean = X{1}(:,indx); %Take particle with highest weight as initial mean

% Map the other particles to a tangent space at wmean
for ii = 1:Npart
    XSampled(:,ii) = logx(xprod(invx(wmean),X{1}(:,ii)));
end

%Compute the mean in the tangent space, check if it is at the origin
Wmean = XSampled*wK';
if norm(Wmean) > 0.01 %If not at the origin, an update is executed
    wmean = xprod(wmean,expx(Wmean)); %New mean on Lie group
end
%Obtain the pose of the weighted mean as output
Y(:,1)= wmean(1:2,:);
Z(:,1)= wmean(1:2,:);


%% 2: UNSCENTED PARTICLE FILTER
for t=2:maxt%length(FileList)
    img=imread(fullfile(append('static/',fname,'/Test_data'), FileList(t).name));                   %Load the image
    hsi = rgb2hsi(img,12,12,4);                                            %Convert image to HSI bins (12x12x4)
    
    %CREATE THE PROPOSAL DISTRIBUTION
    for i=1:Npart %For each particle
        %Run the UKF to create a proposal distribution and the transition prior
        [xPred(:,i),PxxPred{1,i},xPredS{1,i},XPredS{1,i},sigmaN{1,i},weights{1,i},nsp{1,i}]=ukf(X{t-1}(:,i),P{t-1}{:,i},Qv,Rv,alpha,beta,kappa,"CVModel",const); 
        
        %Compute the likelihood of each particle
        [L(i),~,~,~]=likelihood(xPred(:,i),Href,hsi,K,box);
    end    
    
    %OBTAIN A MEASUREMENT
    %Pick the particle with the highest likelihood
    L=L./sum(L); 
    [~,indx]=max(L); %Pick the index
    zmean = xPred(1:2,indx); %Pick the pose (R,o) of the particle with highest likelihood
    for ii = 1:Npart
        % Map the other particles to a tangent space at ZL
        Zpoints(:,ii) = logH(Hprod(invH(zmean),xPred(1:2,ii)));
    end
    
    %Compute the mean in the tangent space, check if it is at the origin
    Zmean = Zpoints*L';     
    while norm(Zmean) > 1e-5 %If not at the origin, an update is executed
        zmean = Hprod(zmean,expH(Zmean)); %New mean on Lie group
        for ii = 1:Npart
            % Map the other particles to a tangent space at wmean
            Zpoints(:,ii) = logH(Hprod(invH(zmean),xPred(1:2,ii)));
        end
        Zmean = Zpoints*L';
    end    
    Z(:,t)= zmean; %Obtain the pose of the weighted mean as output    
    
    %UPDATE THE STATE
    for i=1:Npart
        [xEst(:,i),PEst{1,i}] = updateX(xPred(:,i),PxxPred{1,i},xPredS{1,i},XPredS{1,i},sigmaN{1,i},weights{1,i},nsp{1,i},zmean);
        
        %Sample from the proposal distribution
        xSampled(:,i) = xprod(xEst(:,i),expx(sqrtm(PEst{1,i})*randn(12,1))); 
        
        %Evaluate the sampled particle at the proposal distribution:
        proposal = LieProb(xSampled(:,i),xEst(:,i),PEst{1,i});
        
        %Evaluate the sampled particle at the transition prior:
        prior = LieProb(xSampled(:,i),xPred(:,i),PxxPred{1,i});
        
        %Evaluate the sampled particle at the likelihood function:
        [lik,~,~,~]=likelihood(xSampled(:,i),Href,hsi,K,box);
        
        %Compute the weight of the sampled particle:
        wK(1,i) = lik*prior/proposal;
    end 
    
    %% COMPUTE WEIGHTED MEAN
    wK=wK./sum(wK);     %Normalize the weights
    [~,indx] = max(wK); %Index of the particle with highest weight    
    wmean = xSampled(:,indx); %Take particle with highest weight as initial mean
    
    % Map the other particles to a tangent space at wmean
    for ii = 1:Npart
        XSampled(:,ii) = logx(xprod(invx(wmean),xSampled(:,ii)));
    end
    
    %Compute the mean in the tangent space, check if it is at the origin
    Wmean = XSampled*wK';     
    while norm(Wmean) > 1e-5 %If not at the origin, an update is executed
        wmean = xprod(wmean,expx(Wmean)); %New mean on Lie group
        for ii = 1:Npart
            % Map the other particles to a tangent space at wmean
            XSampled(:,ii) = logx(xprod(invx(wmean),xSampled(:,ii)));
        end
        Wmean = XSampled*wK';
    end
    Y(:,t)= wmean(1:2,:); %Obtain the pose of the weighted mean as output
    
    %% RESAMPLE
    [X{t},count] = Resample_systematic(xSampled,wK);
    P{t} = PEst(1,count);
    %Print progress in command window
    textwaitbar(t, maxt,"Progress   ")
end

%% 3: POST PROCESSING
%Print runtime to command window
runtime = toc(time);
runstring = sprintf('Runtime [s]: %d', runtime);
disp(runstring);

%Save the data
if DoSave
    GUPF_CV_Y = Y;
    save(sprintf(append('Results/',fname,'/GUPF_CV_P%dF%dN1'),Npart,maxt),'GUPF_CV_Y','runtime');
end

%Compute the errors
for ii = 1:maxt
    XEY(ii)  = Y{2,ii}(1)-GT{ii}(1,4); %X-error of the output
    YEY(ii)  = Y{2,ii}(2)-GT{ii}(2,4); %Y-error of the output
    ZEY(ii)  = Y{2,ii}(3)-GT{ii}(3,4); %Z-error of the output
    
    EY = GT{ii}(1:3,4)-Y{2,ii};   
    
    NEY(ii) = norm(EY,1);  %Norm of the error of the output 
    NYR(ii) = rad2deg(norm(logm((GT{ii}(1:3,1:3))\Y{1,ii}))); %||log(R_GT^-1 * R_Y)||%
end 

%% 1: PLOT FIGURES
figure('pos',[200 400 250 200]);
    g2 = plot(XEY,'linewidth',1,'color',[0.9290 0.6940 0.1250]);hold on; grid on;
    xlabel('Frame [-]');
    ylabel('x-error [m]')
    axis([1 maxt -0.15 0.15]);
    XTick = [1 impacts(1) impacts(2) impacts(3) impacts(4) impacts(5) 65];
    XTickLabel = ({'1';num2str(impacts(1));num2str(impacts(2));num2str(impacts(3));num2str(impacts(4));num2str(impacts(5));'65'});
    for ii = 1:length(impacts)
    xline(impacts(ii),':','linewidth',1.2,'color',[0 0 0 1],'alpha',1);
    end
    legend('GUPF\_CV','location','southwest');

figure('pos',[455 400 250 200]);
    plot(YEY,'linewidth',1.0,'color',[0.9290 0.6940 0.1250]);hold on; grid on;
    xlabel('Frame [-]');
    ylabel('y-error [m]')
    axis([1 maxt -0.15 0.15]);
    XTick = [1 impacts(1) impacts(2) impacts(3) impacts(4) impacts(5) 65];
    XTickLabel = ({'1';num2str(impacts(1));num2str(impacts(2));num2str(impacts(3));num2str(impacts(4));num2str(impacts(5));'65'});
    for ii = 1:length(impacts)
    xline(impacts(ii),':','linewidth',1.2,'color',[0 0 0 1],'alpha',1);
    end
    legend('GUPF\_CV','location','southwest');

figure('pos',[710 400 250 200]); 
    plot(ZEY,'linewidth',1.0,'color',[0.9290 0.6940 0.1250]);hold on; grid on;
    xlabel('Frame [-]');
    ylabel('z-error [m]')
    axis([1 maxt -0.15 0.15]);
    XTick = [1 impacts(1) impacts(2) impacts(3) impacts(4) impacts(5) 65];
    XTickLabel = ({'1';num2str(impacts(1));num2str(impacts(2));num2str(impacts(3));num2str(impacts(4));num2str(impacts(5));'65'});
    for ii = 1:length(impacts)
    xline(impacts(ii),':','linewidth',1.2,'color',[0 0 0 1],'alpha',1);
    end
    legend('GUPF\_CV','location','southwest');

figure('pos',[965 400 250 200]); 
    plot(NEY,'linewidth',1.0,'color',[0.9290 0.6940 0.1250]);hold on; grid on;
    xlabel('Frame [-]');
    ylabel('$\|e_{\mathbf{o}}\|$ [m]')
    axis([1 maxt -0.1 0.2]);
    XTick = [1 impacts(1) impacts(2) impacts(3) impacts(4) impacts(5) 65];
    XTickLabel = ({'1';num2str(impacts(1));num2str(impacts(2));num2str(impacts(3));num2str(impacts(4));num2str(impacts(5));'65'});
    for ii = 1:length(impacts)
    xline(impacts(ii),':','linewidth',1.2,'color',[0 0 0 1],'alpha',1);
    end
    legend('GUPF\_CV','location','southwest');

figure('pos',[1220 400 250 200]); 
    plot(NYR,'linewidth',1.0,'color',[0.9290 0.6940 0.1250]);hold on; grid on;
    xlabel('Frame [-]');
    ylabel('$\|e_{\mathbf{R}}\|$ [deg]')
    axis([1 maxt -20 40]);
    XTick = [1 impacts(1) impacts(2) impacts(3) impacts(4) impacts(5) 65];
    XTickLabel = ({'1';num2str(impacts(1));num2str(impacts(2));num2str(impacts(3));num2str(impacts(4));num2str(impacts(5));'65'});
    for ii = 1:length(impacts)
    xline(impacts(ii),':','linewidth',1.2,'color',[0 0 0 1],'alpha',1);
    end
    legend('GUPF\_CV','location','southwest');