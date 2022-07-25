function [channelGain,channelGain_NLoS,HMean_Withoutphase] = RandomAP_generateSetup_Rician_Multi_Antenna(M,K,L,N,nbrOfSetups,correlatedShadowing)

%This function is used to generate realizations of the simulation setup with correlated Rayleigh fading and multi-antenna APs and UEs.
%This is version 1.1 (Last edited: 2022-07-24)

%INPUT:
%M                    = Number of APs for the Cell-free system
%K                    = Number of UEs in the network
%L                    = Number of antennas per AP
%N                    = Number of antennas per UE
%nbrOfSetups          = Number of setups with random UE locations
%correlatedShadowing  = If set to 1 correlated shadow fading is generated, 
%                        else uncorrelated shadow fading
%OUTPUT:
%HMean                = Matrix with dimension MN x K x nbrOfSetups,where (mn',k,n) is the
%                       channel mean between the n'th antenna of AP m and UE k, normalized by
%                       noise power and without random phase shifts
%channelGain          = Matrix with dimension M x K x nbrOfSetups where
%                       element (m,k,n) is the channel gain (normalized by
%                       the noise power) between AP m and UE k in setup n

%% Define simulation setup

%Size of the coverage area (as a square with wrap-around)
cellRange = 1000; %meter

%Communication bandwidth
B = 20e6;

%Noise figure (in dB)
noiseFigure = 7;

%Compute noise power
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Pathloss parameters
alpha_LoS = 26;
constantTerm_LoS = 30.18;

alpha_NLoS = 38;
constantTerm_NLoS = 34.53;

%Standard deviation of the shadow fading
sigma_sf = 8;

%Decorrelation distance of the shadow fading
decorr = 100;

%Shadow fading parameter
delta = 0.5;

%Height difference between an AP and a UE
APheigth = 12.5;  
UEheigth = 1.5;  
VerticalDistance = APheigth-UEheigth;

%The minimum allowed distance (Access Point to User Equipment)
dmin = 10; 

%Define the antenna spacing (in number of wavelengths)
 antennaSpacing = 1/2; %Half wavelength distance
% 
% %Angular standard deviation around the nominal angle (measured in degrees)
% ASDdeg_AP = 15;
% ASDdeg_UE = 15;

%Dropping all UEs while minimum distance requirement is satisfied.
droppedUEs=0;

%Deploy APs randomly
APpositions = cellRange*(rand(M,1) + 1i*rand(M,1));

%Prepare to save results
channelGain_LoS = zeros(M,K,nbrOfSetups);
channelGain_NLoS = zeros(M,K,nbrOfSetups);
HMean_Withoutphase = zeros(M*L,nbrOfSetups,K*N);
%% Go through all setups
for n = 1:nbrOfSetups
 
%------Depoly UEs randomly

%--Compute alternative AP locations by using wrap around
wrapHorizontal = repmat([-cellRange 0 cellRange],[3 1]);
wrapVertical = wrapHorizontal';
wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
APpositionsWrapped = repmat(APpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[M 1]); 

%--Prepare to save the distances and UE positions
distanceAPtoUE = zeros(M,K);
distanceAPtoAP = zeros(M,M); 
distanceUEtoUE = zeros(K,K); 
UEpositions = zeros(K,1);

%--Dropping all UEs
while droppedUEs <K
    
    UEposition = cellRange*(rand(1,1) + 1i*rand(1,1));
    HorizontalDistance = abs(APpositions-UEposition);
    distance = sqrt(VerticalDistance.^2 + HorizontalDistance.^2);
    
    if isempty(distance(distance<dmin))
        droppedUEs = droppedUEs+1;
        distanceAPtoUE(:,droppedUEs)=distance;
       
        %Store UE positions
        UEpositions(droppedUEs)=UEposition;
    end
end


%---Calculate the distances between all AP and UE pairs

%Distances between APs
for m=1:M
    distanceAPtoAP(:,m) = abs(APpositions-APpositions(m));
end

%Distances between UEs
 for k=1:K
    distanceUEtoUE(:,k) = abs(UEpositions-UEpositions(k));
 end
 
%-----Calculate Channel Coefficients

%Create covarince functions for each pair
covMatrixAP = 2.^(-distanceAPtoAP/decorr);
covMatrixUE = 2.^(-distanceUEtoUE/decorr);

%Create shadow fading realizations
shadowFadingAP = sqrtm(covMatrixAP)*randn(M,1); 
shadowFadingUE = sqrt(covMatrixUE)*randn(K,1);

%Create the resulting shadow fading matrix
 shadowFadingMatrix=zeros(M,K);
 for k=1:K
     for m=1:M
        shadowFadingMatrix(m,k)=sqrt(delta)*shadowFadingAP(m) +sqrt(1-delta)*shadowFadingUE(k);
     end
 end
 
 %Scale with variance 
 if correlatedShadowing == 1
    %Correlated Shadow Fading Matrix in dB
    shadowFading = sigma_sf*shadowFadingMatrix; 
 else
    %Uncorrelated Shadow Fading Matrix in dB
    shadowFading = sigma_sf*randn(M,K); 
 end

 
%Prepare to save the result
RicianFactor = zeros(M,K);
channelGaindB = zeros(M,K);
channelGain = zeros(M,K);
%Go through all UEs
 for k = 1:K
     
     [distances_Hori,whichpos] = min(abs(APpositionsWrapped - repmat(UEpositions(k),size(APpositionsWrapped))),[],2);
     distances = sqrt(VerticalDistance^2+distances_Hori.^2);
     
     %Path-loss calculation "Cost 231 Walfish-Ikegami Model"
     %before adding shadow fading
     betaLoS = constantTerm_LoS + alpha_LoS*log10(distances);
     betaNLoS = constantTerm_NLoS + alpha_NLoS*log10(distances);
     
     %Each pair has a LoS
    probLoS = ones(size(distances));%have los
    %   probLoS = zeros(size(distances));%no los
     
     %Calculate the distance based Rician Factor 
      RicianFactor(:,k) = 10.^(1.3-0.003*distances);
     
     %Save the channel gains (in this setup each pair has a LoS path)
     channelGaindB(probLoS==1,k)=-betaLoS(probLoS==1);
     channelGaindB(probLoS==0,k)=-betaNLoS(probLoS==0);
     channelGaindB(:,k) = channelGaindB(:,k)+shadowFading(:,k)-noiseVariancedBm+30;
     channelGain(:,k) = db2pow(channelGaindB(:,k));

     %Scale with Rician factor
     channelGain_LoS(probLoS==1,k,n) = sqrt(RicianFactor(probLoS==1,k)./(RicianFactor(probLoS==1,k) +1 )).*sqrt(channelGain(probLoS==1,k));
     channelGain_NLoS(probLoS==1,k,n) = (1./(RicianFactor(probLoS==1,k) +1 )).*(channelGain(probLoS==1,k));
     channelGain_NLoS(probLoS==0,k,n) = channelGain(probLoS==0,k); %note that probLoS is always one in the manuscript

      for m = 1:M
            
            %Compute nominal angle between UE k and AP m
            angletoUE = angle(UEpositions(k)-APpositionsWrapped(m,whichpos(m)));
            angletoAP = angle(APpositionsWrapped(m,whichpos(m))-UEpositions(k)); 
            ar = zeros(L,1);
            at = zeros(N,1);
            ar(1:L) = exp(1i*2*pi.*(0:(L-1))*sin(angletoUE)*antennaSpacing);
            at(1:N) = exp(1i*2*pi.*(0:(N-1))*sin(angletoAP)*antennaSpacing);

            %Generate HMean_Withoutphase in Multi-Antenna case
            HMean_Withoutphase((m-1)*L+1:m*L,n,(k-1)*N+1:k*N) = channelGain_LoS(m,k)*ar(1:L)*at(1:N)';
          
      end
end
 
end