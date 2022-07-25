function [Hhat_MMSE,F_Pre] = functionChannelEstimates_MMSE(HH_Vec,R_Vec,HMean,F_precoding,nbrOfRealizations,M,K,L,N,tau_p,Pset)

%MMSE channel estimator for Cell-Free setup. The estimation is locally
%performed at the APs.
%Each AP is equipped with L antennas and each UE is equipped with N antennas .
%This is version 1.1 (Last edited: 2022-07-24)         


%Prepare to store MMSE channel estimates
HHhat_MMSE_vec = zeros(L*N,nbrOfRealizations,M,K);
Hhat_MMSE = zeros(M*L,nbrOfRealizations,K*N);
    
%Store identity matrix of size N x N
eyeLN = eye(L*N);

%Generate realizations of normalized noise 
Np = sqrt(0.5)*(randn(N*L,nbrOfRealizations,M,K) + 1i*randn(N*L,nbrOfRealizations,M,K));

%Generate the precoding matrix
F_Pre = zeros(L*N,L*N,K);

for k = 1:K
    
    F_Pre(:,:,k) = kron(F_precoding((k-1)*N+1:k*N,(k-1)*N+1:k*N).',eye(L)); 
    
end
    
for m = 1:M
    for k = 1:K
        
        yp = zeros(L*N,nbrOfRealizations);
        PsiInv = zeros(L*N,L*N);
        inds = Pset(:,k); 

        
        for z = 1:length(inds)
            
            yp = yp + tau_p*F_Pre(:,:,k)*HH_Vec(:,:,m,inds(z));
            PsiInv = PsiInv + tau_p*F_Pre(:,:,inds(z))*R_Vec(:,:,m,inds(z))*F_Pre(:,:,inds(z))';
            
        end
        
        yp = yp + sqrt(tau_p)*Np(:,:,m,k);
        PsiInv = PsiInv + eyeLN;
 
      
        for z = 1:length(inds)
            
            RPsi = R_Vec(:,:,m,inds(z))*F_Pre(:,:,inds(z))'/PsiInv;
            HHhat_MMSE_vec(:,:,m,inds(z)) = RPsi*yp;
            
        end
    end
end


for m = 1:M
    for k = 1:K
        for n = 1:N
            
            Hhat_MMSE((m-1)*L+1:m*L,:,(k-1)*N+n) = HHhat_MMSE_vec((n-1)*L+1:n*L,:,m,k)+HMean((m-1)*L+1:m*L,:,(k-1)*N+n);
            
        end
    end
end
        




