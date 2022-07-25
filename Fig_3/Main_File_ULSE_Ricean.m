warning('off');
clear all
close all
M = 10;
K = 10;
L = [1,2,3,4,5,6,7,8];
N = 4;
tau_c = 200;
nbrOfSetups = 50;
nbrOfRealizations = 750;
%Pilot Reuse factor
w = 1;
%Create the power vector for all UEs (The uplink power is the same(p)at each UE)
px = 0.2;%200 mW
p = px*ones(1,K);
 %Ricean
 SE_Monte_MR_Combining_Level4_2 = zeros(K,nbrOfSetups,length(L));
 SE_Monte_MR_Combining_Level4_wo_2 = zeros(K,nbrOfSetups,length(L));
 
for i = 1:nbrOfSetups
    [channelGain,channelGain_NLoS,channelGain_LoS,UEpositions,APpositions] = RandomAP_generateSetup_Rician_Multi_Antenna_1(M,K,1,1);
   for n = 1:length(L)
   tau_p = K*N/w;   
      [HMean_Withoutphase] = RandomAP_generateSetup_Rician_Multi_Antenna_2(channelGain_LoS,UEpositions,APpositions,M,K,L(n),N,1);
        %with precoding
        [H,HH,HMean,R_Vec,HH_Vec,Omega,F_precoding] = functionChannelGeneration(channelGain_NLoS,HMean_Withoutphase,M,K,N,L(n),nbrOfRealizations,p);
        [Pset] = functionPilotAllocation( R_Vec,M,K,L(n)*N,tau_p/N);
        [Hhat_MMSE,F_Pre] = functionChannelEstimates_MMSE(HH_Vec,R_Vec,HMean,F_precoding,nbrOfRealizations,M,K,L(n),N,tau_p,Pset);
        [C_MMSE_MMSE_Combining] = functionMatrixGeneration(R_Vec,F_Pre,F_precoding,M,K,L(n),N,tau_p,Pset);
        [V_MMSE_Combining] = functionCompute_MMSE_Combining_Matrix(Hhat_MMSE,C_MMSE_MMSE_Combining,nbrOfRealizations,L(n),N,K,M,F_precoding);
        [SE_level4,~] = functionComputeSE_Fully_Centralized_Small_Cell(V_MMSE_Combining,Hhat_MMSE,C_MMSE_MMSE_Combining,F_precoding,tau_c,tau_p,nbrOfRealizations,N,L(n),K,M,1);
        %without precoding
        for k = 1:K
        F_precoding_wo((k-1)*N+1:k*N,(k-1)*N+1:k*N) = sqrt(p(k))*1/sqrt(N)*eye(N);
        end
        [H_wo,HH_wo,HMean_wo,R_Vec_wo,HH_Vec_wo,Omega_wo,~] = functionChannelGeneration(channelGain_NLoS,HMean_Withoutphase,M,K,N,L(n),nbrOfRealizations,p);
        [Pset_wo] = functionPilotAllocation( R_Vec_wo,M,K,L(n)*N,tau_p/N);
        [Hhat_MMSE_wo,F_Pre_wo] = functionChannelEstimates_MMSE(HH_Vec_wo,R_Vec_wo,HMean_wo,F_precoding_wo,nbrOfRealizations,M,K,L(n),N,tau_p,Pset_wo);
        [C_MMSE_MMSE_Combining_wo] = functionMatrixGeneration(R_Vec_wo,F_Pre_wo,F_precoding_wo,M,K,L(n),N,tau_p,Pset_wo);
        [V_MMSE_Combining_wo] = functionCompute_MMSE_Combining_Matrix(Hhat_MMSE_wo,C_MMSE_MMSE_Combining_wo,nbrOfRealizations,L(n),N,K,M,F_precoding_wo);
        [SE_level4_wo,~] = functionComputeSE_Fully_Centralized_Small_Cell(V_MMSE_Combining_wo,Hhat_MMSE_wo,C_MMSE_MMSE_Combining_wo,F_precoding_wo,tau_c,tau_p,nbrOfRealizations,N,L(n),K,M,1);
        disp(['SE with MR combining vector of L4  ' num2str(i)]);
        %Ricean
         SE_Monte_MR_Combining_Level4_2(:,i,n) = SE_level4;
         SE_Monte_MR_Combining_Level4_wo_2(:,i,n) = SE_level4_wo;

    end
    disp([num2str(i) ' AP-Antenna out of ' num2str(L)]);
end

save 'Rician.mat' SE_Monte_MR_Combining_Level4_2 SE_Monte_MR_Combining_Level4_wo_2;