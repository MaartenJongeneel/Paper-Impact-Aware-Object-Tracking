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
DoSave   = false;    %Decide if you want to save 

if config ==1
    %Constant parameters
    impacts = [16 31 39 49 58]; %Frame at which an impact occurs
    
    %Process noise covariance and measurement noise covariance
    Qv = 1e-6*diag([5 5 50 1 1 1 10 10 10 100 100 500]); %[mm]  Process noise covariance: tells me how much i trust my motion model
    Rv = 1e-4*diag([10 10 10 1 1 1]);  %[mm]  Measurement noise covariance: same as the noisy measurement

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
    
    %Process noise covariance and measurement noise covariance
    Qv = 1e-6*diag([5 5 5 1 1 1 10 10 10 100 100 100]); %[mm]  Process noise covariance
    Rv = 1e-4*diag([10 10 10 1 1 1]);                   %[mm]  Measurement noise covariance:

    %Initial state covariance 
    PR = 5e-4*diag([1 1 1]);    Po = 1e-4*diag([1 1 1]);
    Pv = 1e-7*diag([1 1 1]);    Pw = 1e-7*diag([1 1 1]);

    %Initial pose and velocity
    AR_B = GT{1}(1:3,1:3);   %Initial Rotation
    Ao_B = GT{1}(1:3,4);     %Initial position         [m]
    Av_B = [  1; -1.2;   0]; %Initial linear velocity  [m/s]
    Aw_B = [  5;    2;   0]; %Initial angular velocity [rad/s]
end

%Constant parameters
theta = 90;
m     = 1;                  %Mass of the cuboid                [kg]
g     = 9.81;               %Gravitational acceleration        [m/s^2]
te    = 0.0167;             %End time of the simulation (1/60 sec) [s]
Kz_K  = [0;0;1];            %z-component of the K frame
Ky_K  = [0;1;0];            %y-component of the K frame
Kx_K  = [1;0;0];            %x-component of the K frame
l     = box.dim(1);         %Length of the box
w     = box.dim(2);         %Width of the box
h     = box.dim(3);         %Height of the box
const.eN    = 0.2;          %Coefficient of restitution in normal direction       [-]
const.eT    = 0.2;          %Coefficient of restitution in tangential direction   [-]
const.mu    = 0.8;          %Coefficient of friction
const.N     = 100;          %Number of discretization points   [-]
const.dt    = te/const.N;   %Time between two time steps       [s]
const.a     = 0.01;         %Prox point parameter              [-]
const.tol   = 1e-5;         %Error tol for fixed-point         [-]

%Initial guess for LambdaN and LambdaT
LambdaNfull = zeros(8,1);        
LambdaTfull(1:8,1) = {zeros(2,1)}; 

% Force acting on body B: expressed in B with orientation of A:
BA_fo  = [0; m*g; 0];
BA_Tau = [0; 0; 0];
const.BA_f   = [BA_fo; BA_Tau];

%Mass matrix of the cuboid
Ml = m*eye(3);

%Inertia matrix of the cuboid
I = [(m/12)*(w^2+h^2),  0, 0;
     0,  (m/12)*(l^2+h^2), 0;
     0,  0, (m/12)*(l^2+w^2);];

%Inertia tensor
const.B_M_B = [Ml zeros(3,3); zeros(3,3) I];

%Coordinate frame of the contact surface
AR_K    = [1 0 0; 0 cos(deg2rad(theta)) -sin(deg2rad(theta)); 0 sin(deg2rad(theta)) cos(deg2rad(theta))];
const.Ao_K    = [0; 0.4; 1];

const.Az_K = AR_K*Kz_K; %Direction of the normal of the plane in terms of frame A
const.Ay_K = AR_K*Ky_K;
const.Ax_K = AR_K*Kx_K;

%Vertices of the box
const.vertices  = [l -l -l l l -l -l l; -w -w w w -w -w w w; -h -h -h -h h h h h]/2;

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
        [xPred(:,i),PxxPred{1,i},xPredS{1,i},XPredS{1,i},sigmaN{1,i},weights{1,i},nsp{1,i}]=ukf(X{t-1}(:,i),P{t-1}{:,i},Qv,Rv,alpha,beta,kappa,"MotionModel",const); 
        
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
    GUPF_NS_Y = Y;
    save(sprintf(append('Results/',fname,'/GUPF_NS_P%dF%dN1'),Npart,maxt),'GUPF_NS_Y','runtime');
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
    legend('GUPF\_NS','location','southwest');

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
    legend('GUPF\_NS','location','southwest');

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
    legend('GUPF\_NS','location','southwest');

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
    legend('GUPF\_NS','location','southwest');

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
    legend('GUPF\_NS','location','southwest');