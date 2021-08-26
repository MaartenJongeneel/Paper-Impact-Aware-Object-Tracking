close all; clearvars;
addpath(genpath('static')); addpath('Functions'); 
set(groot,'defaulttextinterpreter','latex'); set(groot,'defaultAxesTickLabelInterpreter','latex');  set(groot,'defaultLegendInterpreter','latex');
%% Settings
dosave = 0;
createvideo = 0;
impacts = [14 20 28 35 45]; %Frame at which an impact occurs
ws     = 1;                  %Width of the contact surface   [m]
ls     = 1;                  %Length of the contact surface  [m]
%% Plot results from simulation
FileList = dir(fullfile('GUPF-NS/Test_data', '*.png')); %Load the data set
load('K.mat');       %Camera intrinsic matrix
load('GT.mat');      %Ground truth
load('box.mat');

PF_CV = load('Results/PF_CV_P500F65_5.mat');
PF_NS = load('Results/PF_NS_P500F65N4.mat');
GUPF_NS = load('Results/GUPF_NS_P500F65N3.mat');
% GUPFCV = load('GUPF-CV/Results/50P65FV5.mat');

%% Constants
Kz_K  = [0;0;1];            %z-component of the K frame
Ky_K  = [0;1;0];            %y-component of the K frame
Kx_K  = [1;0;0];            %x-component of the K frame
s     = 1;                  %Size of the bouncing suraface     [m]
theta = 90;                 %Rotation of the bouncing surface [deg]
AR_K  = [1 0 0; 0 cos(deg2rad(theta)) -sin(deg2rad(theta)); 0 sin(deg2rad(theta)) cos(deg2rad(theta))];
Ao_K  = [0; 0.4; 1];
surfacepoints = [0.5*ws -0.5*ws -0.5*ws 0.5*ws 0.5*ws; -0.5*ls -0.5*ls 0.5*ls 0.5*ls -0.5*ls; 0 0 0 0 0;];
spoints = AR_K*surfacepoints +Ao_K;
maxt = 65;
plotYcuboid = 1;

l     = box.dim(1);         %Length of the box
w     = box.dim(2);         %Width of the box
h     = box.dim(3);         %Height of the box
%% Origin of the other frames wrt frame B
Bo_C = [ l/2; -w/2; -h/2];      Bo_G = [ l/2; -w/2;  h/2];
Bo_D = [-l/2; -w/2; -h/2];      Bo_H = [-l/2; -w/2;  h/2];
Bo_E = [-l/2;  w/2; -h/2];      Bo_I = [-l/2;  w/2;  h/2];
Bo_F = [ l/2;  w/2; -h/2];      Bo_J = [ l/2;  w/2;  h/2];

%% Loading results
Y{1} = PF_CV.PF_CV_Y;
Y{2} = PF_NS.PF_NS_Y;
Y{3} = GUPF_NS.GUPF_NS_Y;

for t = 1:maxt
q{1}(:,t)  = Y{1}{2,t};   %Position of Y
q{2}(:,t)  = Y{2}{2,t};   %Position of Y
q{3}(:,t)  = Y{3}{2,t};   %Position of Y

R1{1}(:,t) = Y{1}{1,t}(1:3,1);   %Rotation around x of Y
R2{1}(:,t) = Y{1}{1,t}(1:3,2);   %Rotation around y of Y
R3{1}(:,t) = Y{1}{1,t}(1:3,3);   %Rotation around z of Y

R1{2}(:,t) = Y{2}{1,t}(1:3,1);   %Rotation around x of Y
R2{2}(:,t) = Y{2}{1,t}(1:3,2);   %Rotation around y of Y
R3{2}(:,t) = Y{2}{1,t}(1:3,3);   %Rotation around z of Y

R1{3}(:,t) = Y{3}{1,t}(1:3,1);   %Rotation around x of Y
R2{3}(:,t) = Y{3}{1,t}(1:3,2);   %Rotation around y of Y
R3{3}(:,t) = Y{3}{1,t}(1:3,3);   %Rotation around z of Y

end
for t = 1:65
GTo(:,t)  = GT{t}(1:3,4);         %Postion of GT
GTR1(:,t) = GT{t}(1:3,1);        %Rotation around x of GT
GTR2(:,t) = GT{t}(1:3,2);        %Rotation around y of GT
GTR3(:,t) = GT{t}(1:3,3);        %Rotation around z of GT
end 



for ii = 1:maxt
    %Compute the errors of each Y w.r.t. the ground truth
    EY{1}(:,ii) = GT{ii}(1:3,4)-q{1}(:,ii); 
    EY{2}(:,ii) = GT{ii}(1:3,4)-q{2}(:,ii);
    EY{3}(:,ii) = GT{ii}(1:3,4)-q{3}(:,ii);
    
    NEY{1}(ii) = norm(EY{1}(:,ii),1);  %Norm of the error of the output
    NEY{2}(ii) = norm(EY{2}(:,ii),1);  %Norm of the error of the output
    NEY{3}(ii) = norm(EY{3}(:,ii),1);  %Norm of the error of the output
    
    NYR{1}(ii) = rad2deg(norm(logm((GT{ii}(1:3,1:3))\Y{1}{1,ii}))); %||log(R_GT^-1 * R_Z)||%
    NYR{2}(ii) = rad2deg(norm(logm((GT{ii}(1:3,1:3))\Y{2}{1,ii}))); %||log(R_GT^-1 * R_Z)||%
    NYR{3}(ii) = rad2deg(norm(logm((GT{ii}(1:3,1:3))\Y{3}{1,ii}))); %||log(R_GT^-1 * R_Z)||%
end

%% Plot figures
close all 

ple = [1 2 3];
p1 = AR_K*[ 1/2*2*s;  1/2*2*s; 0];
p2 = AR_K*[ 1/2*2*s; -1/2*2*s; 0];
p3 = AR_K*[-1/2*2*s; -1/2*2*s; 0];
p4 = AR_K*[-1/2*2*s;  1/2*2*s; 0];

xco = [p1(1) p2(1) p3(1) p4(1)]+Ao_K(1);
yco = [p1(3) p2(3) p3(3) p4(3)]+Ao_K(3);
zco = [p1(2) p2(2) p3(2) p4(2)]+Ao_K(2);

%Plot figures at these xy coordinates
px = [10 265 520 810 1065 1355 1540];
py = [30 320 610 900]+15;
for  ii = 1:length(px)
    for jj = 1:length(py)
        pp{jj,ii} = [px(ii) py(jj)];
    end 
end 

RMS_PF_CV   = rms(NEY{ple(1)});
RMS_PF_NS = rms(NEY{ple(2)});
RMS_GUPF_NS = rms(NEY{ple(3)});

EredPF_NS = (1-(RMS_PF_NS/RMS_PF_CV))*100;
EredGUPF_NS = (1-(RMS_GUPF_NS/RMS_PF_CV))*100;

R_RMS_PF_CV   = rms(NYR{ple(1)});
R_RMS_PF_NS = rms(NYR{ple(2)});
R_RMS_GUPF_NS = rms(NYR{ple(3)});

R_EredPF_NS = (1-(R_RMS_PF_NS/R_RMS_PF_CV))*100;
R_EredGUPF_NS = (1-(R_RMS_GUPF_NS/R_RMS_PF_CV))*100;

%% ------------------------------- FIGURE 1 ------------------------------- %%

traj = 3; %choose the trajectory to plot
figure('rend','painters','pos',[pp{2,2} 250 200]); 
    ha = tight_subplot(1,1,[.08 .07],[.1 .03],[0.1 0.0]); %[gap_h gap_w] [low up ] [lft rght]
    axes(ha(1));
    
    %Plot the camera
    cam = plotCamera('Location',[0 0 0],'Orientation',[1 0 0; 0 0 -1; 0 1 0],'Opacity',0,'Size',0.05,'color',[0.2 0.2 0.2]);
    
    count = 1;
    frameplot = [1 6 12 18 24 30 36 42 48 54 60 65];
    plot3(GTo(1,:),GTo(3,:),GTo(2,:),'color','k');  hold on;   %Plot GT trajectory
    plot3(q{traj}(1,:),q{traj}(3,:),q{traj}(2,:),'color','b'); %Plot Y trajectory
 
    for ii=1:maxt
        if frameplot(count) == ii
        
        %Create the cuboid
        AR_B = Y{traj}{1,ii};
        Ao_B = Y{traj}{2,ii};
        Ao_C = AR_B*Bo_C + Ao_B;
        Ao_D = AR_B*Bo_D + Ao_B;
        Ao_E = AR_B*Bo_E + Ao_B;
        Ao_F = AR_B*Bo_F + Ao_B;
        Ao_G = AR_B*Bo_G + Ao_B;
        Ao_H = AR_B*Bo_H + Ao_B;
        Ao_I = AR_B*Bo_I + Ao_B;
        Ao_J = AR_B*Bo_J + Ao_B;
        
        %Plot the cuboid of Y
        plot3([Ao_C(1) Ao_D(1)],[Ao_C(3) Ao_D(3)],[Ao_C(2) Ao_D(2)],'b');%
        plot3([Ao_D(1) Ao_E(1)],[Ao_D(3) Ao_E(3)],[Ao_D(2) Ao_E(2)],'b');%
        plot3([Ao_E(1) Ao_F(1)],[Ao_E(3) Ao_F(3)],[Ao_E(2) Ao_F(2)],'b');
        plot3([Ao_F(1) Ao_C(1)],[Ao_F(3) Ao_C(3)],[Ao_F(2) Ao_C(2)],'b');
        plot3([Ao_G(1) Ao_H(1)],[Ao_G(3) Ao_H(3)],[Ao_G(2) Ao_H(2)],'b');%
        plot3([Ao_H(1) Ao_I(1)],[Ao_H(3) Ao_I(3)],[Ao_H(2) Ao_I(2)],'b');%
        plot3([Ao_I(1) Ao_J(1)],[Ao_I(3) Ao_J(3)],[Ao_I(2) Ao_J(2)],'b');
        plot3([Ao_J(1) Ao_G(1)],[Ao_J(3) Ao_G(3)],[Ao_J(2) Ao_G(2)],'b');
        plot3([Ao_C(1) Ao_G(1)],[Ao_C(3) Ao_G(3)],[Ao_C(2) Ao_G(2)],'b');
        plot3([Ao_D(1) Ao_H(1)],[Ao_D(3) Ao_H(3)],[Ao_D(2) Ao_H(2)],'b');
        plot3([Ao_E(1) Ao_I(1)],[Ao_E(3) Ao_I(3)],[Ao_E(2) Ao_I(2)],'b');
        plot3([Ao_F(1) Ao_J(1)],[Ao_F(3) Ao_J(3)],[Ao_F(2) Ao_J(2)],'b');
        
        %Create the cuboid of GT
        AR_B = GT{ii}(1:3,1:3);
        Ao_B = GT{ii}(1:3,4);
        Ao_C = AR_B*Bo_C + Ao_B;
        Ao_D = AR_B*Bo_D + Ao_B;
        Ao_E = AR_B*Bo_E + Ao_B;
        Ao_F = AR_B*Bo_F + Ao_B;
        Ao_G = AR_B*Bo_G + Ao_B;
        Ao_H = AR_B*Bo_H + Ao_B;
        Ao_I = AR_B*Bo_I + Ao_B;
        Ao_J = AR_B*Bo_J + Ao_B;
        
        %Plot the cuboid of Y
        plot3([Ao_C(1) Ao_D(1)],[Ao_C(3) Ao_D(3)],[Ao_C(2) Ao_D(2)],'k');%
        plot3([Ao_D(1) Ao_E(1)],[Ao_D(3) Ao_E(3)],[Ao_D(2) Ao_E(2)],'k');%
        plot3([Ao_E(1) Ao_F(1)],[Ao_E(3) Ao_F(3)],[Ao_E(2) Ao_F(2)],'k');
        plot3([Ao_F(1) Ao_C(1)],[Ao_F(3) Ao_C(3)],[Ao_F(2) Ao_C(2)],'k');
        plot3([Ao_G(1) Ao_H(1)],[Ao_G(3) Ao_H(3)],[Ao_G(2) Ao_H(2)],'k');%
        plot3([Ao_H(1) Ao_I(1)],[Ao_H(3) Ao_I(3)],[Ao_H(2) Ao_I(2)],'k');%
        plot3([Ao_I(1) Ao_J(1)],[Ao_I(3) Ao_J(3)],[Ao_I(2) Ao_J(2)],'k');
        plot3([Ao_J(1) Ao_G(1)],[Ao_J(3) Ao_G(3)],[Ao_J(2) Ao_G(2)],'k');
        plot3([Ao_C(1) Ao_G(1)],[Ao_C(3) Ao_G(3)],[Ao_C(2) Ao_G(2)],'k');
        plot3([Ao_D(1) Ao_H(1)],[Ao_D(3) Ao_H(3)],[Ao_D(2) Ao_H(2)],'k');
        plot3([Ao_E(1) Ao_I(1)],[Ao_E(3) Ao_I(3)],[Ao_E(2) Ao_I(2)],'k');
        plot3([Ao_F(1) Ao_J(1)],[Ao_F(3) Ao_J(3)],[Ao_F(2) Ao_J(2)],'k');
        
        count = count+1;
        end
    end
    %Plot the inclined table K
    table3 = fill3(spoints(1,1:4),spoints(3,1:4),spoints(2,1:4),1);hold on;
    set(table3,'FaceColor',[0.8 0.8 0.8],'FaceAlpha',1);

    grid on;axis equal;%axis off;
    axis([-0.6 0.6 -0.05 1.7 -0.1 0.4]);
    xlabel('x [m]','position',[  0 -0.2 0.5]);
    ylabel('z [m]','position',[  -0.7 0.88 0.5]);
    zlabel('y [m]','position',[  -0.86 1.6 0]);
    set(gca, 'ZDir','reverse')
    view(-26,34);
    hold off
    legend('GT','GUPF\_NS','location','northeast');
    if dosave ==1
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,'Figures/GUPF_NS_Trajectory.pdf','-dpdf','-painters')
    end


%% ------------------------------- FIGURE 2 ------------------------------- %%
figure('rend','painters','pos',[pp{3,2} 1.2*450 200]); 
    ha = tight_subplot(1,3,[.08 .09],[.18 .18],[0.08 0.02]);  %[gap_h gap_w] [lower upper] [left right]
    axes(ha(1));
    plot(EY{ple(1)}(1,:),'-','linewidth',1,'color',[0.9290 0.6940 0.1250]);hold on; grid on;
    plot(EY{ple(2)}(1,:),'-.','linewidth',1,'color','b');
    plot(EY{ple(3)}(1,:),'-.','linewidth',1.5,'color','r');
    xlabel('Frame [-]');
    ylabel('$x$-error, $e_x$ [m]')
    axis([1 maxt -0.1 0.1]);
    ha(1).XTick = [1 impacts(1) impacts(2) impacts(3) impacts(4) impacts(5) 65];
    ha(1).XTickLabel = ({'1';num2str(impacts(1));num2str(impacts(2));num2str(impacts(3));num2str(impacts(4));num2str(impacts(5));'65'});
    for ii = 1:length(impacts)
    xline(impacts(ii),':','linewidth',1.2,'color',[0 0 0 1],'alpha',1);
    end
    
    axes(ha(2));
    plot(EY{ple(1)}(2,:),'-','linewidth',1,'color',[0.9290 0.6940 0.1250]);hold on; grid on;
    plot(EY{ple(2)}(2,:),'-.','linewidth',1,'color','b');
    plot(EY{ple(3)}(2,:),'-.','linewidth',1.5,'color','r');
    xlabel('Frame [-]');
    ylabel('$y$-error, $e_y$ [m]')
    set(gca, 'YDir','reverse')
    axis([1 maxt -0.1 0.1]);
    ha(2).XTick = [1 impacts(1) impacts(2) impacts(3) impacts(4) impacts(5) 65];
    ha(2).XTickLabel = ({'1';num2str(impacts(1));num2str(impacts(2));num2str(impacts(3));num2str(impacts(4));num2str(impacts(5));'65'});
    for ii = 1:length(impacts)
    xline(impacts(ii),':','linewidth',1.2,'color',[0 0 0 1],'alpha',1);
    end
   
    axes(ha(3));
    plot(EY{ple(1)}(3,:),'-','linewidth',1,'color',[0.9290 0.6940 0.1250]);hold on; grid on;
    plot(EY{ple(2)}(3,:),'-.','linewidth',1,'color','b');
    plot(EY{ple(3)}(3,:),'-.','linewidth',1.5,'color','r');
    xlabel('Frame [-]');
    ylabel('$z$-error, $e_z$ [m]')
    axis([1 maxt -0.10 0.1]);   
    ha(3).XTick = [1 impacts(1) impacts(2) impacts(3) impacts(4) impacts(5) 65];
    ha(3).XTickLabel = ({'1';num2str(impacts(1));num2str(impacts(2));num2str(impacts(3));num2str(impacts(4));num2str(impacts(5));'65'});
    for ii = 1:length(impacts)
    xline(impacts(ii),':','linewidth',1.2,'color',[0 0 0 1],'alpha',1);
    end
    
    L1 = legend('PF\_CV','PF\_NS','GUPF\_NS','Impact Time','NumColumns',4,'location','northeast');
    L1.Position(2) = 0.88;
    L1.Position(1) = 0.5-(L1.Position(3)/2);
    if dosave ==1
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,'Figures/POSErrors.pdf','-dpdf','-painters')
    end   
    
    
%% ------------------------------- FIGURE 3 ------------------------------- %%

figure('rend','painters','pos',[pp{2,3} 1.2*450 200]);
    ha = tight_subplot(1,2,[0 .08],[.18 .18],[0.09 0.02]);%[gap_h gap_w] [low up ] [lft rght]
    
    axes(ha(1));
    plot(NEY{ple(1)},'-','linewidth',1,'color',[0.9290 0.6940 0.1250]);hold on; grid on;
    plot(NEY{ple(2)},'-.','linewidth',1,'color','b');
    plot(NEY{ple(3)},'-.','linewidth',1.5,'color','r');
    plot([1 65],[0 0],'-','linewidth',1,'color','k');
    xlabel('Frame [-]');
    ylabel('Position error $\|e_{\mathbf{o}}\|$ [m]')
    axis([1 maxt -0.05 0.25]);    
    ha(1).XTick = [1 impacts(1) impacts(2) impacts(3) impacts(4) impacts(5) 65];
    ha(1).XTickLabel = ({'1';num2str(impacts(1));num2str(impacts(2));num2str(impacts(3));num2str(impacts(4));num2str(impacts(5));'65'});
    for ii = 1:length(impacts)
    xline(impacts(ii),':','linewidth',1.2,'color',[0 0 0 1],'alpha',1);
    end
    
    axes(ha(2));
    g1= plot(NYR{ple(1)},'-','linewidth',1,'color',[0.9290 0.6940 0.1250]);hold on; grid on;
    g2= plot(NYR{ple(2)},'-.','linewidth',1,'color','b');
    g3= plot(NYR{ple(3)},'-.','linewidth',1.5,'color','r');
    plot([1 65],[0 0],'-','linewidth',1,'color','k');
    xlabel('Frame [-]');
    ylabel('Rotation error $\|e_{\mathbf{R}}\|$ [deg]')
    axis([1 maxt -4 20]);
    ha(2).XTick = [1 impacts(1) impacts(2) impacts(3) impacts(4) impacts(5) 65];
    ha(2).XTickLabel = ({'1';num2str(impacts(1));num2str(impacts(2));num2str(impacts(3));num2str(impacts(4));num2str(impacts(5));'65'});
    for ii = 1:length(impacts)
    g4 = xline(impacts(ii),':','linewidth',1.2,'color',[0 0 0 1],'alpha',1);
    end
    
    L1 = legend([g1 g2 g3 g4],'PF\_CV','PF\_NS','GUPF\_NS','Impact Time','NumColumns',4,'location','northeast');
    L1.Position(2) = 0.88;
    L1.Position(1) = 0.5-(L1.Position(3)/2);
    if dosave ==1
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,'Figures/POS_ROT_Errors.pdf','-dpdf','-painters')
    end
    

%% ------------------------------- FIGURE 4 ------------------------------- %%

figure('rend','painters','pos',[pp{3,4} 1.2*450 200]); 
    ha = tight_subplot(1,3,[.08 .09],[.18 .18],[0.08 0.02]);  %[gap_h gap_w] [lower upper] [left right]
    axes(ha(1));
    plot(GTo(1,:),'-','linewidth',2,'color','k');hold on; grid on;
    plot(q{ple(1)}(1,:),'-','linewidth',1,'color',[0.9290 0.6940 0.1250]);
    plot(q{ple(2)}(1,:),'-.','linewidth',1,'color','b');
    plot(q{ple(3)}(1,:),'-.','linewidth',1.5,'color','r');
    xlabel('Frame [-]');
    ylabel('Position x [m]')
    axis([1 maxt -0.3 0.4]);
    ha(1).XTick = [1 impacts(1) impacts(2) impacts(3) impacts(4) impacts(5) 65];
    ha(1).XTickLabel = ({'1';num2str(impacts(1));num2str(impacts(2));num2str(impacts(3));num2str(impacts(4));num2str(impacts(5));'65'});
    for ii = 1:length(impacts)
    xline(impacts(ii),':','linewidth',1.2,'color',[0 0 0 1],'alpha',1);
    end
    
    axes(ha(2));
    plot(GTo(2,:),'-','linewidth',2,'color','k');hold on; grid on;
    plot(q{ple(1)}(2,:),'-','linewidth',1,'color',[0.9290 0.6940 0.1250]);
    plot(q{ple(2)}(2,:),'-.','linewidth',1,'color','b');
    plot(q{ple(3)}(2,:),'-.','linewidth',1.5,'color','r');
    xlabel('Frame [-]');
    ylabel('Position y [m]')
    set(gca, 'YDir','reverse')
    axis([1 maxt -0.05 0.45]);
    ha(2).XTick = [1 impacts(1) impacts(2) impacts(3) impacts(4) impacts(5) 65];
    ha(2).XTickLabel = ({'1';num2str(impacts(1));num2str(impacts(2));num2str(impacts(3));num2str(impacts(4));num2str(impacts(5));'65'});
    for ii = 1:length(impacts)
    xline(impacts(ii),':','linewidth',1.2,'color',[0 0 0 1],'alpha',1);
    end
   
    axes(ha(3));
    plot(GTo(3,:),'-','linewidth',2,'color','k');hold on; grid on;
    plot(q{ple(1)}(3,:),'-','linewidth',1,'color',[0.9290 0.6940 0.1250]);
    plot(q{ple(2)}(3,:),'-.','linewidth',1,'color','b');
    plot(q{ple(3)}(3,:),'-.','linewidth',1.5,'color','r');
    xlabel('Frame [-]');
    ylabel('Position z [m]')
    axis([1 maxt 0.75 1.52]);
    ha(3).XTick = [1 impacts(1) impacts(2) impacts(3) impacts(4) impacts(5) 65];
    ha(3).XTickLabel = ({'1';num2str(impacts(1));num2str(impacts(2));num2str(impacts(3));num2str(impacts(4));num2str(impacts(5));'65'});
    for ii = 1:length(impacts)
    xline(impacts(ii),':','linewidth',1.2,'color',[0 0 0 1],'alpha',1);
    end   
    
    L1 = legend('GT','PF\_CV','PF\_NS','GUPF\_NS','Impact Time','NumColumns',5,'location','northeast');
    L1.Position(2) = 0.88;
    L1.Position(1) = 0.5-(L1.Position(3)/2);
    if dosave ==1
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,'Figures/ResultingTrajectory.pdf','-dpdf','-painters')
    end
   
    
%% ------------------------------- FIGURE 5 ------------------------------- %%

figure('rend','painters','pos',[pp{3,6} 450 200]); 
    ha = tight_subplot(1,3,[.08 -0.2],[-0.05 0],[-0.1 -0.1]);  %[gap_h gap_w] [lower upper] [left right]
    axes(ha(1));    
    [Xding,Yding,Zding] = sphere(20);   % draw the pi-ball
    hSurface = surf(Xding,Yding,Zding);hold on;
    set(hSurface,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.3,'FaceLighting','gouraud','EdgeColor','none','LineStyle','none','BackFaceLighting','unlit',...
                 'AmbientStrengt',0.3,'DiffuseStrength',0.7,'SpecularStrength',0.5,'SpecularExponent',5,'SpecularColorReflectance',1);
    plot3(GTR1(1,1),GTR1(2,1),GTR1(3,1),'*','MarkerSize',10,'color','k');
    
    plot3(GTR1(1,1),GTR1(2,1),GTR1(3,1),'*','MarkerSize',10,'color','k');
    plot3(GTR1(1,1:49),GTR1(2,1:49),GTR1(3,1:49),'color','k','linewidth',2.5);
    plot3(GTR1(1,49:65),GTR1(2,49:65),GTR1(3,49:65),'color',[0 0 0 0.1],'linewidth',2.5);

    plot3(R1{ple(1)}(1,1:49),R1{ple(1)}(2,1:49),R1{ple(1)}(3,1:49),'color',[0.9290 0.6940 0.1250],'linewidth',1.2);
    plot3(R1{ple(1)}(1,49:65),R1{ple(1)}(2,49:65),R1{ple(1)}(3,49:65),'color',[0.9290 0.6940 0.1250 0.1],'linewidth',1.2);
    
    plot3(R1{ple(2)}(1,1:49),R1{ple(2)}(2,1:49),R1{ple(2)}(3,1:49),'--','color','b','linewidth',1.2);
    plot3(R1{ple(2)}(1,49:65),R1{ple(2)}(2,49:65),R1{ple(2)}(3,49:65),'--','color',[0 0 1 0.1],'linewidth',1.2);
    
    plot3(R1{ple(3)}(1,1:49),R1{ple(3)}(2,1:49),R1{ple(3)}(3,1:49),'--','color','r','linewidth',1.2);
    plot3(R1{ple(3)}(1,49:65),R1{ple(3)}(2,49:65),R1{ple(3)}(3,49:65),'--','color',[1 0 0 0.1],'linewidth',1.2);
    
    view(128,31)
    light('Position',[-1 1 1])
    axis equal
    axis off
    text(0,0,-1.5,'x')

    axes(ha(2));    
    [Xding,Yding,Zding] = sphere(20);   % draw the pi-ball
    hSurface = surf(Xding,Yding,Zding);hold on;
    set(hSurface,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.2,'FaceLighting','gouraud','EdgeColor','none','LineStyle','none','BackFaceLighting','unlit',...
                 'AmbientStrengt',0.3,'DiffuseStrength',0.7,'SpecularStrength',0.5,'SpecularExponent',5,'SpecularColorReflectance',1);
    plot3(GTR2(1,1),GTR2(2,1),GTR2(3,1),'*','MarkerSize',7,'color','k');
    plot3(GTR2(1,1:7),GTR2(2,1:7),GTR2(3,1:7),'color','k','linewidth',2.5);
    plot3(GTR2(1,7:20),GTR2(2,7:20),GTR2(3,7:20),'color',[0 0 0 0.3],'linewidth',2.5);
    plot3(GTR2(1,20:65),GTR2(2,20:65),GTR2(3,20:65),'color','k','linewidth',2.5);

    plot3(R2{ple(1)}(1,1:7),R2{ple(1)}(2,1:7),R2{ple(1)}(3,1:7),'color',[0.9290 0.6940 0.1250],'linewidth',1.2);
    plot3(R2{ple(1)}(1,7:20),R2{ple(1)}(2,7:20),R2{ple(1)}(3,7:20),'color',[0.9290 0.6940 0.1250 0.5],'linewidth',1.2);
    plot3(R2{ple(1)}(1,20:65),R2{ple(1)}(2,20:65),R2{ple(1)}(3,20:65),'color',[0.9290 0.6940 0.1250],'linewidth',1.2);
    
    plot3(R2{ple(2)}(1,1:7),R2{ple(2)}(2,1:7),R2{ple(2)}(3,1:7),'--','color','b','linewidth',1.2);
    plot3(R2{ple(2)}(1,7:20),R2{ple(2)}(2,7:20),R2{ple(2)}(3,7:20),'--','color',[0 0 1 0.5],'linewidth',1.2);
    plot3(R2{ple(2)}(1,20:65),R2{ple(2)}(2,20:65),R2{ple(2)}(3,20:65),'--','color','b','linewidth',1.2);
    
    plot3(R2{ple(3)}(1,1:7),R2{ple(3)}(2,1:7),R2{ple(3)}(3,1:7),'--','color','r','linewidth',1.2);
    plot3(R2{ple(3)}(1,7:20),R2{ple(3)}(2,7:20),R2{ple(3)}(3,7:20),'--','color',[1 0 0 0.5],'linewidth',1.2);
    plot3(R2{ple(3)}(1,20:65),R2{ple(3)}(2,20:65),R2{ple(3)}(3,20:65),'--','color','r','linewidth',1.2);
    
    view(128,31)
    light('Position',[-1 1 1])
    axis equal
    axis off
    text(0,0,-1.5,'y')

    axes(ha(3));    
    [Xding,Yding,Zding] = sphere(20);   % draw the pi-ball
    hSurface = surf(Xding,Yding,Zding);hold on;
    set(hSurface,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.2,'FaceLighting','gouraud','EdgeColor','none','LineStyle','none','BackFaceLighting','unlit',...
                 'AmbientStrengt',0.3,'DiffuseStrength',0.7,'SpecularStrength',0.5,'SpecularExponent',5,'SpecularColorReflectance',1);
    plot3(GTR3(1,1),GTR3(2,1),GTR3(3,1),'*','MarkerSize',7,'color','k');
    plot3(GTR3(1,1:15),GTR3(2,1:15),GTR3(3,1:15),'color',[0 0 0 0.3],'linewidth',2.5);
    g1 = plot3(GTR3(1,15:65),GTR3(2,15:65),GTR3(3,15:65),'color','k','linewidth',2.5);

    plot3(R3{ple(1)}(1,1:15),R3{ple(1)}(2,1:15),R3{ple(1)}(3,1:15),'color',[0.9290 0.6940 0.1250 0.5],'linewidth',1.2);
    g2 = plot3(R3{ple(1)}(1,15:65),R3{ple(1)}(2,15:65),R3{ple(1)}(3,15:65),'color',[0.9290 0.6940 0.1250],'linewidth',1.2);
    
    plot3(R3{ple(2)}(1,1:15),R3{ple(2)}(2,1:15),R3{ple(2)}(3,1:15),'--','color',[0 0 1 0.5],'linewidth',1.2);
    g3 =plot3(R3{ple(2)}(1,15:65),R3{ple(2)}(2,15:65),R3{ple(2)}(3,15:65),'--','color','b','linewidth',1.2);
    
    plot3(R3{ple(3)}(1,1:15),R3{ple(3)}(2,1:15),R3{ple(3)}(3,1:15),'--','color',[1 0 0 0.5],'linewidth',1.2);
    g4 =plot3(R3{ple(3)}(1,15:65),R3{ple(3)}(2,15:65),R3{ple(3)}(3,15:65),'--','color','r','linewidth',1.2);
    
    view(128,31)
    light('Position',[-1 1 1])
    axis equal
    axis off
    text(0,0,-1.5,'z')

    L1 = legend([g1 g2 g3 g4],{'GT','PF\_CV','PF\_NS','GUPF\_NS'},'NumColumns',4,'location','northeast');
    L1.Position(2) = 0.88;
    L1.Position(1) = 0.5-(L1.Position(3)/2);
    L1.FontSize = 6;
    
    if dosave ==1
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,'Figures/ResultsRotations.pdf','-dpdf','-painters')
    end

    
%% ------------------------------- FIGURE 6 ------------------------------- %%  
    
figure('rend','painters','pos',[pp{2,5} 400 200]); 
    ha = tight_subplot(1,1,[.08 .07],[.18 .1],[0.12 0.02]);  %[gap_h gap_w] [lower upper] [left right]
    axes(ha(1));
    plot(GTo(2,:),'-','linewidth',2,'color','k');hold on; grid on;
    plot(q{ple(1)}(2,:),'-','linewidth',1,'color',[0.9290 0.6940 0.1250]);
    plot(q{ple(2)}(2,:),'-.','linewidth',1,'color','b');
    plot(q{ple(3)}(2,:),'-.','linewidth',1.5,'color','r');
    xlabel('Frame [-]');
    ylabel('y [m]')
    set(gca, 'YDir','reverse')
    axis([impacts(1)-3 impacts(2)+3 0.23 0.36]);    
    ha(1).XTick = [1 impacts(1) impacts(2) impacts(3) impacts(4) impacts(5) 65];
    ha(1).XTickLabel = ({'1';num2str(impacts(1));num2str(impacts(2));num2str(impacts(3));num2str(impacts(4));num2str(impacts(5));'65'});
    for ii = 1:length(impacts)
    xline(impacts(ii),':','linewidth',1.2,'color',[0 0 0 1],'alpha',1);
    end
    L1 = legend('GT','PF\_CV','PF\_NS','GUPF\_NS','Impact Time','location','eastoutside');
    if dosave ==1
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,'Figures/CloseUpY.pdf','-dpdf','-painters')
    end
    
    
  