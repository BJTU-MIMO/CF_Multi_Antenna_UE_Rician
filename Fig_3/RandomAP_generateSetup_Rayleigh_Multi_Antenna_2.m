function [HMean_Withoutphase] = RandomAP_generateSetup_Rayleigh_Multi_Antenna_2(channelGain_LoS,UEpositions,APpositions,M,K,L,N,nbrOfSetups)
cellRange = 1000; 
wrapHorizontal = repmat([-cellRange 0 cellRange],[3 1]);
wrapVertical = wrapHorizontal';
wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
APpositionsWrapped = repmat(APpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[M 1]); 
antennaSpacing = 0.5;
HMean_Withoutphase = zeros(M*L,nbrOfSetups,K*N);

for k = 1:K 
    [distances_Hori,whichpos] = min(abs(APpositionsWrapped - repmat(UEpositions(k),size(APpositionsWrapped))),[],2);
for m = 1:M
            
            %Compute nominal angle between UE k and AP m
            angletoUE = angle(UEpositions(k)-APpositionsWrapped(m,whichpos(m)));
            angletoAP = angle(APpositionsWrapped(m,whichpos(m))-UEpositions(k)); 
            ar = zeros(L,1);
            at = zeros(N,1);
            ar(1:L) = exp(1i*2*pi.*(0:(L-1))*sin(angletoUE)*antennaSpacing);
            at(1:N) = exp(1i*2*pi.*(0:(N-1))*sin(angletoAP)*antennaSpacing);

            %Generate HMean_Withoutphase in Multi-Antenna case
            HMean_Withoutphase((m-1)*L+1:m*L,1,(k-1)*N+1:k*N) = channelGain_LoS(m,k)*ar(1:L)*at(1:N)';
          
      end
end