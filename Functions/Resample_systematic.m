function [out,count] = Resample_systematic(in,wn)
% Systematic resampling of the particles in by the weights wn
% INPUTS:    in           : Particles: double(4,Npart) column is a particle
%            wn           : Normalized weights of the particles
%
% OUTPUTS:   out          : Resampled particles
%% SYSTEMATIC RESAMPLING:
[r,N] = size(in);                %N = Number of particles.
out = cell(r,N);                 %Preallocate memory
count = zeros(1,N);              %Preallocate memory
Q=cumsum(wn);                    %Commulative sum of weights
T=linspace(0,1-1/N,N)+rand(1,N)./N; %A number is chosen according to:
                                    %uk_s = ((k-1)+u*)/N with u*~U[0,1)
i=1;                    %Count in T
j=1;                    %Count in Q
while i<=N
    while(T(i)>Q(j))    %As long as T (current spoke) is bigger than Q(current weight)
        j=j+1;          %go to the next Q value
    end
    out(:,i) = in(:,j); %When Q (the weight) is bigger, copy the particle to the output
    count(i)=j;
    i=i+1;              %Go to next T(spoke) value  
end


