warning('off');
clear all
close all
load Rayleigh.mat;
load Rician.mat;

K=10;
nbrOfSetups=50;
L=[1,2,3,4,5,6,7,8];
M1 = zeros(length(L),1);
M2 = zeros(length(L),1);
M3 = zeros(length(L),1);
M4 = zeros(length(L),1);
for n = 1:length(L)
    MMM1 = SE_Monte_MR_Combining_Level4_1(:,:,n);
    MMM2 = SE_Monte_MR_Combining_Level4_wo_1(:,:,n);
    MM1 = sort(MMM1(:));
    MM2 = sort(MMM2(:));
    M1(n) = MM1(K*nbrOfSetups*0.1);
    M2(n) = MM2(K*nbrOfSetups*0.1);
    MMM3 = SE_Monte_MR_Combining_Level4_2(:,:,n);
    MMM4 = SE_Monte_MR_Combining_Level4_wo_2(:,:,n);
    MM3 = sort(MMM3(:));
    MM4 = sort(MMM4(:));
    M3(n) = MM3(K*nbrOfSetups*0.1);
    M4(n) = MM4(K*nbrOfSetups*0.1);
end
c1 = plot(L,M1,'r-s','LineWidth',1.3);
hold on
c2 = plot(L,M2,'k-o','LineWidth',1.3,'MarkerSize',5);
hold on
c3 = plot(L,M3,'r-s','LineWidth',1.3);
hold on
c4 = plot(L,M4,'k-o','LineWidth',1.3,'MarkerSize',5);
grid on
grid minor
xlabel('Number of AP antennas, L','Interpreter','Latex');
ylabel('95%-likely per-user UL SE [bit/s/Hz]','Interpreter','Latex');
set(gca,'FontSize',12);
legend([c1 c2],{'precoding','without precoding'},'Interpreter','Latex','Location','Northwest');