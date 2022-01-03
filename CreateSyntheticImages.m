%% Create synthetic images of the box for a given trajectory
% This script creates the synthetic images. It uses the BoxSimulator
% function to simulate the trajectory of the box, given an (initial) input 
% pose and velocity of the box. Under the settings below, one can determine
% this inital pose and velocity, together with some characteristics of the
% box, like dimensions, mass, and coefficients used in the simulator, as
% well as defining the contact surface and the camera intrinsics. One can
% choose to store the data in one of the config folders under /static, or
% create a new one. The synthetic images will be stored there, together
% with the camera intrinsic matrix, ground truth data of the box, and the
% box model that were used to create the images. 

%% General stuff
close all; clearvars; clc; addpath('Functions'); addpath('images');
%% Settings
releaseOrientation = Rx(20);        %Release orientation of the box            [deg]
releasePosition    = [0.5; 0.8; 0]; %Release position of the box               [m]
releaseLinVel      = [-2; 0; 0];    %Release linear velocity (expressed in B)  [m/s]
releaseAngVel      = [5; 2; 0];     %Release angular velocity (expressed in B) [rad/s]
eN                 = 0.2;           %Normal coefficient of restitution         [-]
eT                 = 0.2;           %Tangential coefficient of restitution     [-]
mu                 = 0.8;           %Coefficient of friction                   [-]
l                  = 0.1;           %Length of the box                         [m]
w                  = 0.15;          %Width of the box                          [m]
h                  = 0.05;          %Height of the box                         [m]
m                  = 1;             %Mass of the box                           [kg]
AR_C               = Rx(0);         %Orientation of the contact plane          [deg]
Ao_C               = [0; 1; -0.4];  %Position of the contact plane             [m]
s                  = 0.6;           %size of the contact plane (square)        [m]
runtime            = 1.1;           %Runtime of the simulation                 [s]
fps                = 1/60;          %Framerate for the video output            [fps]
dt                 = fps/10;        %Timestep at which the simulator runs      [s]
doPlot             = true;          %Decide if we want to plot the box         [-]
createvideo        = true;          %Decide if we want to create a video       [-]
configFolder       = 'static/config03'; %Config folder where images are stored [-]

%Camera intrinsic matrix
K = [400    0  320;
    0  400  240;
    0    0    1];

%Load the background image
img = imread('grey.jpg');

%% Do some checks
if isfolder(configFolder)
    str = input('Configuration folder already exists, do you want to overwrite the data? [y/n]','s');
    if contains(str,'y')
        disp(append('Data of config folder ',configFolder,' will be overwritten.'));
        dosave = true;
    elseif contains(str,'n')
        disp('Data will not be stored');
        dosave = false;
    else
        error(append('Invalid input given. You entered: ',str));
    end
else
    str = input('Configuration folder does not yet exist, do you want to create it and save the data? [y/n]','s');
    if contains(str,'y')
        mkdir(configFolder);
        mkdir(append(configFolder,'/Test_data'));
        disp(append('New config folder ',configFolder,' has been made, data will be stored there'));
        dosave = true;
    elseif contains(str,'n')
        disp('No folder has been created, no data will be stored');
        dosave = false;
    else
        error(append('Invalid input given. You entered: ',str));
    end
end

%% Create the box model and run the simulation
box = create_box_model(l,w,h,m,doPlot);
AH_B = BoxSimulator(releasePosition,releaseOrientation,releaseLinVel,releaseAngVel,eN,eT,mu,box,AR_C,Ao_C,dt,runtime);

%% Plotting and saving the image
tel = 1;
for ii = 1:(fps/dt):length(AH_B)
    %Obtain the state of the object form its pose
    pose = [Rx(90) zeros(3,1); zeros(1,3),1]*AH_B{ii};
    state{1} = pose(1:3,1:3);
    state{2} = pose(1:3,4);
    
    %Show the resulting image and save the image
    imgout = ProjectionCuboid(img,K,box,state,-Ao_C(3),Ao_C(2),s,640,480);
    
    %Plot the results to the screen
    if doPlot
        figure(1);
        imshow(imgout);hold on;
        hold off;
    end
    
    %Save the results to
    if dosave
        imwrite(imgout,fullfile(append(configFolder,'/Test_data/'), sprintf('%03d.png',tel)))
    end
    
    %Get the ground truth pose
    GT{tel} = pose;
    tel=tel+1;
end

if dosave
    save(append(configFolder,'/GT.mat'),'GT');
    save(append(configFolder,'/K.mat'),'K');
    save(append(configFolder,'/box.mat'),'box');
end

%% Create a video of the box
if createvideo
    video = VideoWriter(append(configFolder,'/BouncingVideoBig'),'MPEG-4'); %create the video object
    video.FrameRate = 30;
    open(video); %open the file for writing
    for ii = 1:(fps/dt):length(AH_B)
        %Obtain the state of the object form its pose
        pose = [Rx(90) zeros(3,1); zeros(1,3),1]*AH_B{ii};
        state{1} = pose(1:3,1:3);
        state{2} = pose(1:3,4);
        
        %Show the resulting image and save the image
        imgout = ProjectionCuboid(img,K,box,state,-Ao_C(3),Ao_C(2),s,640,480);
        writeVideo(video,imgout); %write the image to file
    end
    close(video); %close the file
end