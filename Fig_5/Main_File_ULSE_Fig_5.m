warning('off');
clear all
close all
M = [10,20,40,80];
K = 10;
L = 4;
N = [1,2,3,4,5,6];
tau_c = 200;
nbrOfSetups = 50;
nbrOfRealizations = 750;
%Pilot Reuse factor
w = 1;
%Create the power vector for all UEs (The uplink power is the same(p)at each UE)
px = 0.2;%200 mW
p = px*ones(1,K);
 
 SE_Monte_MMSE_Combining_Level4 = zeros(K,nbrOfSetups,length(N),length(M));
for m=1:length(M)
for i = 1:nbrOfSetups
       [channelGain,channelGain_NLoS,channelGain_LoS,UEpositions,APpositions] = RandomAP_generateSetup_Rician_Multi_Antenna_1(M(m),K,1,1); 
      for n = 1:length(N)  
         tau_p = K*N(n)/w; 
         [HMean_Withoutphase] = RandomAP_generateSetup_Rician_Multi_Antenna_2(channelGain_LoS,UEpositions,APpositions,M(m),K,L,N(n),1);

        [H,HH,HMean,R_Vec,HH_Vec,Omega,F_precoding] = functionChannelGeneration(channelGain_NLoS,HMean_Withoutphase,M(m),K,N(n),L,nbrOfRealizations,p);
        [Pset] = functionPilotAllocation( R_Vec,M(m),K,L*N(n),tau_p/N(n));
        [Hhat_MMSE,F_Pre] = functionChannelEstimates_MMSE(HH_Vec,R_Vec,HMean,F_precoding,nbrOfRealizations,M(m),K,L,N(n),tau_p,Pset);
        [C_MMSE_MMSE_Combining] = functionMatrixGeneration(R_Vec,F_Pre,F_precoding,M(m),K,L,N(n),tau_p,Pset);
        [V_MMSE_Combining] = functionCompute_MMSE_Combining_Matrix(Hhat_MMSE,C_MMSE_MMSE_Combining,nbrOfRealizations,L,N(n),K,M(m),F_precoding);
        [SE_level4,~] = functionComputeSE_Fully_Centralized_Small_Cell(V_MMSE_Combining,Hhat_MMSE,C_MMSE_MMSE_Combining,F_precoding,tau_c,tau_p,nbrOfRealizations,N(n),L,K,M(m),0);  

        SE_Monte_MMSE_Combining_Level4(:,i,n,m) = SE_level4;
        disp([num2str(i) ' AP-Antenna out of ' num2str(L)]);
      end
       
end
end

for i=1:length(N)
    MM1 = SE_Monte_MMSE_Combining_Level4(:,:,i,1);
    MM2 = SE_Monte_MMSE_Combining_Level4(:,:,i,2);
    MM3 = SE_Monte_MMSE_Combining_Level4(:,:,i,3);
    MM4 = SE_Monte_MMSE_Combining_Level4(:,:,i,4);
    M1(i) = mean(MM1(:));
    M2(i) = mean(MM2(:));
    M3(i) = mean(MM3(:));
    M4(i) = mean(MM4(:));
end
figure;
c1=plot(N,M1,'k-o','LineWidth',1.3);
hold on;
c2=plot(N,M2,'b-+','LineWidth',1.3);
hold on;
c3=plot(N,M3,'m-s','LineWidth',1.3);
hold on;
c4=plot(N,M4,'r->','LineWidth',1.3);
grid on
grid minor
xlabel('Number of antennas per UE $(N)$','Interpreter','Latex');
ylabel('Average UL SE [bit/s/Hz]','Interpreter','Latex');
set(gca,'FontSize',12);
legend([c4,c3,c2,c1],{'M=80','M=40','M=20','M=10'},'Interpreter','Latex','Location','Northwest');