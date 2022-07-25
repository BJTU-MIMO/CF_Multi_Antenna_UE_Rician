function [H,HH,HMean,R_Vec,HH_Vec,Omega,F_precoding] = functionChannelGeneration(channelGain_NLoS,HMean_Withoutphase,M,K,N,L,nbrOfRealizations,p)

%---This function is used to generate the Weichselberger correlation model. The outputs are channel realizations and channel means with
%random phase shifts at a coherence block.
%Each AP is equipped with L antennas and each UE is equipped with N antennas.
%This is version 1.1 (Last edited: 2022-07-24)

%Prepare to store the results   
U_APR = zeros(M*L,M*L,K);
U_UET = zeros(K*N,K*N,M);

%---UL Precoding
F_precoding = zeros(K*N,K*N);
Omega = zeros(M*L,K*N);

WW = sqrt(0.5)*(randn(M*L,K*N,nbrOfRealizations)+1i*randn(M*L,K*N,nbrOfRealizations));

R_Vec = zeros(L*N,L*N,M,K);
HH_Vec = zeros(L*N,nbrOfRealizations,M,K);
HH = zeros(M*L,nbrOfRealizations,K*N);
HMean = zeros(M*L,nbrOfRealizations,K*N);

for k = 1:K
    for m = 1:M
        
        %----Unitary Matrix
        [U_APR((m-1)*L+1:m*L,(m-1)*L+1:m*L,k),~,~] = svd(rand(L,L) + 1i*rand(L,L));
        [U_UET((k-1)*N+1:k*N,(k-1)*N+1:k*N,m),~,~] = svd(rand(N,N) + 1i*rand(N,N)); 

        %----Coupling Matrix
        c = 0.8 + 0.15*rand(1);
        a1 = rand(L,1);
        a11 = a1./sum(a1(:))*L*N*c;
        a111 = sort(a11,'descend');
        a2 = rand(L,N-1);%---N>1
        a22 = a2./sum(a2(:))*L*N*(1-c);
        Omega0 = [a111,a22];
        
        Omega((m-1)*L+1:m*L,(k-1)*N+1:k*N) = channelGain_NLoS(m,k)*Omega0; %--
           
    end
end


%---Channel Generation
Hp = permute((sqrt(Omega).*WW),[1,3,2]);
for k = 1:K
    for i = (k-1)*N+1:k*N
        
        HH(:,:,i) = U_APR(:,:,k)*Hp(:,:,i);
        
    end
end

for m = 1:M
    for t = (m-1)*L+1:m*L
        
        H11 = reshape(HH(t,:,:),nbrOfRealizations,K*N);
        HH(t,:,:) = reshape(H11*U_UET(:,:,m)',nbrOfRealizations,K*N);
        
    end
end
%---phaseshift
angles = -pi + 2*pi*rand(M,nbrOfRealizations,K);
phaseMatrix = exp(1i*angles);
v_LN = ones(L,N);
PhaseMatrix = zeros(M*L,nbrOfRealizations,K*N);
for i = 1:nbrOfRealizations
    for m = 1:M
        for k = 1:K
            PhaseMatrix((m-1)*L+1:m*L,i,(k-1)*N+1:k*N) = phaseMatrix(m,i,k)*v_LN;
        end
    end
    
    HMean(:,i,:) = PhaseMatrix(:,i,:).*HMean_Withoutphase;
end
 H = HMean + HH;
%---
%---Full Correlation Matrix & Channel Vectorization
HH_v = permute(HH,[1,3,2]);
for m = 1:M
    for k = 1:K
        
        Omeg = reshape(Omega((m-1)*L+1:m*L,(k-1)*N+1:k*N),L*N,1);
        R_Vec(:,:,m,k) = kron(conj(U_UET((k-1)*N+1:k*N,(k-1)*N+1:k*N,m)),U_APR((m-1)*L+1:m*L,(m-1)*L+1:m*L,k))*diag(Omeg)*kron(conj(U_UET((k-1)*N+1:k*N,(k-1)*N+1:k*N,m)),U_APR((m-1)*L+1:m*L,(m-1)*L+1:m*L,k))';
        HH_Vec(:,:,m,k) = reshape(HH_v((m-1)*L+1:m*L,(k-1)*N+1:k*N,:),L*N,nbrOfRealizations);
        
    end
end


UU_t = sum(U_UET,3);

for k = 1:K
    
    F_precoding((k-1)*N+1:k*N,(k-1)*N+1:k*N) = sqrt(p(k))*1/norm(UU_t((k-1)*N+1:k*N,(k-1)*N+1:k*N),'fro')*UU_t((k-1)*N+1:k*N,(k-1)*N+1:k*N);
    
end

