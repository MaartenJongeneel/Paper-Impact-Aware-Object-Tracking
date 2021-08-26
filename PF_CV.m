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
DoSave   = true;    %Decide if you want to save results

if config ==1
    %Constant parameters
    impacts = [16 31 39 49 58]; %Frame at which an impact occurs
    const.Noise = 0.035*diag([1.1 1.1 0.85 0.9 0.9 0.9]);

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
    impacts = [14 20 28 35 45]; %Frame at which an impact occurs
    const.Noise = 0.001*diag([1.7*1 9.9*1.2 1*1 0.7*24*5 0.7*12*2 0.7*6*0.1]);

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
img=imread(fullfile(append('static/',fname,'/Test_data'), FileList(1).name));  %Load the image
hsi = rgb2hsi(img,12,12,4);                           %Convert image to HSI bins (12x12x4)
for i = 1:Npart
    xiR = zeros(3,1) + sqrtm(PR)*randn(3,1);
    xio = zeros(3,1) + sqrtm(Po)*randn(3,1);
    xiv = zeros(3,1) + sqrtm(Pv)*randn(3,1);
    xiw = zeros(3,1) + sqrtm(Pw)*randn(3,1);
    
    X{1}{1,i}  = AR_B*expm(hat(xiR));
    X{1}{2,i}  = Ao_B+xio;
    X{1}{3,i}  = Av_B+xiv;
    X{1}{4,i}  = Aw_B+xiw;
             
    [wK(1,i),~,~,~]=likelihood(X{1}(:,i),Href,hsi,K,box); 
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

%% 2: PARTICLE FILTER
for t=2:maxt
    img=imread(fullfile(append('static/',fname,'/Test_data'), FileList(t).name)); %Load the image
    hsi = rgb2hsi(img,12,12,4);                          %Convert image to HSI bins (12x12x4)
        
    % PROPAGATE THE PARTICLES (CREATE PRIOR/PROPOSAL) AND COMPUTE LIKELIHOOD
    for i=1:Npart %For each particle
        X{t}(:,i) = CVModel(X{t-1}(:,i),const);
        [wK(1,i),~,~,~]=likelihood(X{t}(:,i),Href,hsi,K,box);
    end
    
    % COMPUTE WEIGHTED MEAN
    wK=wK./sum(wK);       %Normalize the weights
    [~,indx] = max(wK);   %Index of the particle with highest weight
    wmean = X{t}(:,indx); %Take particle with highest weight as initial mean
    
    % Map the other particles to a tangent space at wmean
    for ii = 1:Npart
        XSampled(:,ii) = logx(xprod(invx(wmean),X{t}(:,ii)));
    end
    
    %Compute the mean in the tangent space, check if it is at the origin
    Wmean = XSampled*wK'; 
    while norm(Wmean) > 1e-5 %If not at the origin, an update is executed
        wmean = xprod(wmean,expx(Wmean)); %New mean on Lie group
        for ii = 1:Npart
            % Map the other particles to a tangent space at wmean
            XSampled(:,ii) = logx(xprod(invx(wmean),X{t}(:,ii)));
        end
        Wmean = XSampled*wK';
    end
    %Obtain the pose of the weighted mean as output
    Y(:,t)= wmean(1:2,:); 
    
    % RESAMPLE
    [X{t},count] = Resample_systematic(X{t},wK);
    %Print progress in command window
    textwaitbar(t, maxt,"Progress   ")
end

%% 3: POSTPROCESSING
%Print runtime to command window
runtime = toc(time);
runstring = sprintf('Runtime [s]: %d', runtime);
disp(runstring);

%Save the data
if DoSave
    PF_CV_Y = Y;
    save(sprintf(append('Results/',fname,'/PF_CV_P%dF%dN2'),Npart,maxt),'PF_CV_Y','runtime');
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
close all;

% Figures to plot positions and errors
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
    legend('PF\_CV','location','southwest');

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
    legend('PF\_CV','location','southwest');

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
    legend('PF\_CV','location','southwest');

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
    legend('PF\_CV','location','southwest');

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
    legend('PF\_CV','location','southwest');