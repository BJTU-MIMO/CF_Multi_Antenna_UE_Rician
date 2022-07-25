warning('off');
clear all
close all
M = [10,20,30,40,50,60,70,80,90,100];
K = 10;
L = 2;
N = 2;
tau_c = 200;
nbrOfSetups = 50;
nbrOfRealizations = 750;
%Pilot Reuse factor
w = 2;
%Create the power vector for all UEs (The uplink power is the same(p)at each UE)
px = 0.2;%200 mW
p = px*ones(1,K);
 
 SE_Monte_MMSE_Combining_Level3 = zeros(K,nbrOfSetups,length(M));
 SE_Monte_MR_Combining_Level3 = zeros(K,nbrOfSetups,length(M));
 SE_Monte_MMSE_Combining_Level4 = zeros(K,nbrOfSetups,length(M));
 SE_Monte_MR_Combining_Level4 = zeros(K,nbrOfSetups,length(M));
 SE_Th_MR_Combining_Level3 = zeros(K,nbrOfSetups,length(M)); 
 
for m = 1:length(M)
    tau_p = K*N/w;
    for i = 1:nbrOfSetups

         [channelGain,channelGain_NLoS,HMean_Withoutphase] = RandomAP_generateSetup_Rician_Multi_Antenna(M(m),K,L,N,1,1);
         [H,HH,HMean,R_Vec,HH_Vec,Omega,F_precoding] = functionChannelGeneration(channelGain_NLoS,HMean_Withoutphase,M(m),K,N,L,nbrOfRealizations,p);
         [Pset] = functionPilotAllocation( R_Vec,M(m),K,L*N,tau_p/N);
         [Hhat_MMSE,F_Pre] = functionChannelEstimates_MMSE(HH_Vec,R_Vec,HMean,F_precoding,nbrOfRealizations,M(m),K,L,N,tau_p,Pset);

         [SE_MR_Combining_LSFD] = functionComputeMonteCarloSE_UL_LSFD(Hhat_MMSE,HH,HMean,tau_c,tau_p,nbrOfRealizations,N,L,K,M(m),F_precoding);
         disp(['SE with MR combining vector of L3  ' num2str(i)]);
         [C_MMSE_MMSE_Combining] = functionMatrixGeneration(R_Vec,F_Pre,F_precoding,M(m),K,L,N,tau_p,Pset);
         [V_MMSE_Combining] = functionCompute_MMSE_Combining_Matrix(Hhat_MMSE,C_MMSE_MMSE_Combining,nbrOfRealizations,L,N,K,M(m),F_precoding);
         [SE_MMSE_Combining_LSFD] = functionComputeMonteCarloSE_UL_LSFD(V_MMSE_Combining,HH,HMean,tau_c,tau_p,nbrOfRealizations,N,L,K,M(m),F_precoding);
         disp(['SE with MMSE combining vector of L3  ' num2str(i)]);
         [SE_MR_Combining_L4,~] = functionComputeSE_Fully_Centralized_Small_Cell(Hhat_MMSE,Hhat_MMSE,C_MMSE_MMSE_Combining,F_precoding,tau_c,tau_p,nbrOfRealizations,N,L,K,M(m),1);
         disp(['SE with MR combining vector of L4  ' num2str(i)]);
         [SE_MMSE_Combining_L4,~] = functionComputeSE_Fully_Centralized_Small_Cell(V_MMSE_Combining,Hhat_MMSE,C_MMSE_MMSE_Combining,F_precoding,tau_c,tau_p,nbrOfRealizations,N,L,K,M(m),0);
         disp(['SE with MMSE combining vector of L4  ' num2str(i)]);
         [SE_MR_th_LSFD] = functionCompute_Theoretical_SE( R_Vec,F_precoding,F_precoding,HMean_Withoutphase,M(m),K,L,N,tau_p,tau_c,Pset);
         disp(['Theoretical SE with MR combining vector of L3  ' num2str(i)]);

         SE_Monte_MMSE_Combining_Level3(:,i,m) = SE_MMSE_Combining_LSFD;
         SE_Monte_MR_Combining_Level3(:,i,m) = SE_MR_Combining_LSFD;
         SE_Monte_MMSE_Combining_Level4(:,i,m) = SE_MMSE_Combining_L4;
         SE_Monte_MR_Combining_Level4(:,i,m) = SE_MR_Combining_L4;
         SE_Th_MR_Combining_Level3(:,i,m) = SE_MR_th_LSFD; 
         
    end
    disp([num2str(m) ' AP-Antenna out of ' num2str(L)]);
    MM1 = SE_Monte_MMSE_Combining_Level3(:,:,m);
    MM2 = SE_Monte_MR_Combining_Level3(:,:,m);
    MM3 = SE_Monte_MMSE_Combining_Level4(:,:,m);
    MM4 = SE_Monte_MR_Combining_Level4(:,:,m);
    MM5 = SE_Th_MR_Combining_Level3(:,:,m);
    M1(m) = mean(MM1(:));
    M2(m) = mean(MM2(:));
    M3(m) = mean(MM3(:));
    M4(m) = mean(MM4(:));
    M5(m) = mean(MM5(:));
end

figure;
c1=plot(M,M3,'b-s','LineWidth',1.3);
hold on;
c2=plot(M,M4,'b--s','LineWidth',1.3);
hold on;
c3=plot(M,M1,'r-o','LineWidth',1.3);
hold on
c4=plot(M,M2,'r--o','LineWidth',1.3);
hold on
c5=plot(M,M5,'kx','LineWidth',1.3);
hold on;
grid on
grid minor
xlabel('Number of APs $(M)$','Interpreter','Latex');
ylabel('Average UL SE [bit/s/Hz]','Interpreter','Latex');
set(gca,'FontSize',12);
legend([c1,c2,c3,c4,c5],{'FCP (MMSE)','FCP (MR)','LSFD (L-MMSE)','LSFD (MR)','Analytical results'},'Interpreter','Latex','Location','Northwest');