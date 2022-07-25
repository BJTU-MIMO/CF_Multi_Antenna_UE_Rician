warning('off');
clear all
close all
M = 20;
K = 10;
L = [1,2,3,4,5,6,7,8];
N = 2;
tau_c = 200;
nbrOfSetups = 50;
nbrOfRealizations = 750;
%Pilot Reuse factor
w = 1;
%Create the power vector for all UEs (The uplink power is the same(p)at each UE)
px = 0.2;%200 mW
p = px*ones(1,K);
 
 SE_Monte_MR_Combining_Level3 = zeros(K,nbrOfSetups,length(L));
 SE_Monte_MR_Combining_Level4 = zeros(K,nbrOfSetups,length(L));
 SE_Monte_MR_Combining_Level2 = zeros(K,nbrOfSetups,length(L));
 
for i = 1:nbrOfSetups
       [channelGain,channelGain_NLoS,channelGain_LoS,UEpositions,APpositions] = RandomAP_generateSetup_Rician_Multi_Antenna_1(M,K,1,1); 
      for l = 1:length(L)  
         tau_p = K*N/w; 
         [HMean_Withoutphase] = RandomAP_generateSetup_Rician_Multi_Antenna_2(channelGain_LoS,UEpositions,APpositions,M,K,L(l),N,1);

        [H,HH,HMean,R_Vec,HH_Vec,Omega,F_precoding] = functionChannelGeneration(channelGain_NLoS,HMean_Withoutphase,M,K,N,L(l),nbrOfRealizations,p);
        [Pset] = functionPilotAllocation( R_Vec,M,K,L(l)*N,tau_p/N);
        [Hhat_MMSE,F_Pre] = functionChannelEstimates_MMSE(HH_Vec,R_Vec,HMean,F_precoding,nbrOfRealizations,M,K,L(l),N,tau_p,Pset);
         [C_MMSE_MMSE_Combining] = functionMatrixGeneration(R_Vec,F_Pre,F_precoding,M,K,L(l),N,tau_p,Pset);
         [SE_level2] = functionComputeMonteCarloSE_UL_L2(Hhat_MMSE,HH,HMean,tau_c,tau_p,nbrOfRealizations,N,L(l),K,M,F_precoding);
         [SE_level3] = functionComputeMonteCarloSE_UL_L3(Hhat_MMSE,HH,HMean,tau_c,tau_p,nbrOfRealizations,N,L(l),K,M,F_precoding);
         [SE_level4,~] = functionComputeSE_Fully_Centralized_Small_Cell(Hhat_MMSE,Hhat_MMSE,C_MMSE_MMSE_Combining,F_precoding,tau_c,tau_p,nbrOfRealizations,N,L(l),K,M,1);  

         SE_Monte_MR_Combining_Level3(:,i,l) = SE_level3;
         SE_Monte_MR_Combining_Level4(:,i,l) = SE_level4;
         SE_Monte_MR_Combining_Level2(:,i,l) = SE_level2;
        disp([num2str(l) ' AP-Antenna out of ' num2str(L)]);
      end
       
end
for l=1:length(L)
    MM1 = SE_Monte_MR_Combining_Level4(:,:,l);
    MM2 = SE_Monte_MR_Combining_Level3(:,:,l);
    MM3 = SE_Monte_MR_Combining_Level2(:,:,l);
    M1(l) = mean(MM1(:));
    M2(l) = mean(MM2(:));
    M3(l) = mean(MM3(:));
end
figure;
c1=plot(L,M1,'b-s','LineWidth',1.3);
hold on;
c2=plot(L,M2,'r-o','LineWidth',1.3);
hold on;
c3=plot(L,M3,'k->','LineWidth',1.3);
hold on;
grid on
grid minor
xlabel('Number of antennas per AP $(L)$','Interpreter','Latex');
ylabel('Average UL SE [bit/s/Hz]','Interpreter','Latex');
set(gca,'FontSize',12);
legend([c1,c2,c3],{'FCP','LSFD','SCD'},'Interpreter','Latex','Location','Northwest');
