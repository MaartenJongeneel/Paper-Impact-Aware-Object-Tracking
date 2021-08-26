function [X1] = MotionModel(X0,const)
% Computes the state of the cuboid for the next time step
% INPUTS:    X0        : State at t0 cell(R;o;v;w)
%
% OUTPUTS:   X1        : State at t1 cell(R;o;v;w)
%% Compute the state of the next time step
%% Initial pose and velocity
AR_B    = X0{1};
Ao_B    = X0{2};
AH_B = [AR_B, Ao_B; zeros(1,3), 1];

BV_AB = [X0{3}; X0{4}];

LambdaNfull = zeros(8,1);   %Initial guess for LambdaN        
LambdaTfull(1:8,1) = {zeros(2,1)}; %initial guess for LambdaT

%Retrieve constants from struct
Bp_1 = const.vertices(:,1); Bp_5 = const.vertices(:,5);
Bp_2 = const.vertices(:,2); Bp_6 = const.vertices(:,6);
Bp_3 = const.vertices(:,3); Bp_7 = const.vertices(:,7);
Bp_4 = const.vertices(:,4); Bp_8 = const.vertices(:,8);

B_M_B = const.B_M_B;
%% Dynamics
for t = 1:const.N %For each time step
    %Kinematics: Compute the configuration at the mid-time
    AH_Bm = AH_B*expm(0.5*const.dt*hat(BV_AB));
    AR_Bm = AH_Bm(1:3,1:3);
    Ao_Bm = AH_Bm(1:3,4);
    AR_B  = AH_B(1:3,1:3);
    
    %Compute the wrench at the mid-time
    B_fM = ([AR_Bm zeros(3); zeros(3) AR_Bm])'*const.BA_f;
    
    %And compute the gap-functions at the mid-time
    gNmC = const.Az_K'*((Ao_Bm + AR_Bm*Bp_1)-const.Ao_K);
    gNmD = const.Az_K'*((Ao_Bm + AR_Bm*Bp_2)-const.Ao_K);
    gNmE = const.Az_K'*((Ao_Bm + AR_Bm*Bp_3)-const.Ao_K);
    gNmF = const.Az_K'*((Ao_Bm + AR_Bm*Bp_4)-const.Ao_K);
    gNmG = const.Az_K'*((Ao_Bm + AR_Bm*Bp_5)-const.Ao_K);
    gNmH = const.Az_K'*((Ao_Bm + AR_Bm*Bp_6)-const.Ao_K);
    gNmI = const.Az_K'*((Ao_Bm + AR_Bm*Bp_7)-const.Ao_K);
    gNmJ = const.Az_K'*((Ao_Bm + AR_Bm*Bp_8)-const.Ao_K);
    
    gN = [gNmC;gNmD;gNmE;gNmF;gNmG;gNmH;gNmI;gNmJ]; %column of normal contact distances
   
    %Obtain the linear and angular velocity at tA
    vA = BV_AB;
    
    IN = find(gN<0);
    if  IN > 0
        %Compute the matrix containing the normal force directions at ta and tm
        WNta = CompWN(AR_B,const.Az_K,Bp_1,Bp_2,Bp_3,Bp_4,Bp_5,Bp_6,Bp_7,Bp_8);       
        WNtm = CompWN(AR_Bm,const.Az_K,Bp_1,Bp_2,Bp_3,Bp_4,Bp_5,Bp_6,Bp_7,Bp_8); 
        %Compute the matrix containing the tangential force directions at ta and tm
        WTta = CompWT(AR_B,const.Ax_K,const.Ay_K,Bp_1,Bp_2,Bp_3,Bp_4,Bp_5,Bp_6,Bp_7,Bp_8);  
        WTtm = CompWT(AR_Bm,const.Ax_K,const.Ay_K,Bp_1,Bp_2,Bp_3,Bp_4,Bp_5,Bp_6,Bp_7,Bp_8); 
        
        WNA = WNta(:,IN);
        WNM = WNtm(:,IN);
        WTA = cell2mat(WTta(:,IN));
        WTM = cell2mat(WTtm(:,IN));
        converged = 0;
        LambdaN=LambdaNfull(IN);
        LambdaT=cell2mat(LambdaTfull(IN));
        while converged==0
            %Decompose the system to write the linear and angular velocity
            %in different equations
            vE = vA + B_M_B\(B_fM*const.dt - [hat(vA(4:6)), zeros(3); hat(vA(1:3)), hat(vA(4:6))]*B_M_B*vA*const.dt + WNM*LambdaN + WTM*LambdaT);
            
            %Define the normal velocities at the beginning and end of the
            %time step
            gammaNA = WNA'*vA;
            gammaNE = WNM'*vE;
            
            gammaTA = WTA'*vA;
            gammaTE = WTM'*vE;
            
            %Newtons restitution law
            xiN = gammaNE+const.eN*gammaNA;
            xiT = gammaTE+const.eT*gammaTA;
            
            %Find LambdaN using the proximal point function
            LambdaNold = LambdaN;
            LambdaTold = LambdaT;
            LambdaN = proxCN(LambdaN-const.a*xiN);
            LambdaT = proxCT(LambdaT-const.a*xiT,const.mu*LambdaN);
            
            error= norm(LambdaN-LambdaNold)+norm(LambdaT-LambdaTold);
            converged = error<const.tol;
        end
        BV_AB = vE;
    else
        %Update the velocity to the next time step
        vE = B_M_B\(B_fM*const.dt - [hat(vA(4:6)), zeros(3); hat(vA(1:3)), hat(vA(4:6))]*B_M_B*vA*const.dt) + vA;
        BV_AB = vE;
        LambdaN = 0;
        LambdaT = [0;0];
    end
    %Update Lambda for next estimate
    if IN ~= 0
        LambdaNfull(IN)=LambdaN;
        cnt=1;
        for ii = length(IN)
            LambdaTfull(IN(ii)) = {LambdaT(cnt:cnt+1)};
            cnt=cnt+2;
        end
    end
    
    %Complete the time step
    AH_B  = AH_Bm*expm(0.5*const.dt*hat(BV_AB));
end
    
    X1{1,1} = AH_B(1:3,1:3);
    X1{2,1} = AH_B(1:3,4);
    X1{3,1} = BV_AB(1:3);
    X1{4,1} = BV_AB(4:6);
