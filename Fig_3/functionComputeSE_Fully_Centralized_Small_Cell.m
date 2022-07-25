function [SE_level4,SE_l1] = functionComputeSE_Fully_Centralized_Small_Cell(H_Combining,H_hat,C,F_precoding,tau_c,tau_p,nbrOfRealizations,N,L,K,M,V_d)
warning('off');
%Compute uplink SE for Cell-free mMIMO for the four different receiver
%cooperation levels, using either MR or MMSE/L-MMSE combining
%This is version 1.1 (Last edited: 2022-07-24)

%Store identity matrices of different sizes
eyeL = eye(L);
eyeML = eye(M*L);

%Compute the prelog factor
prelogFactor = (1-tau_p/tau_c);

%Compute sum of all estimation error correlation matrices at every BS
C_tot = zeros(L,L,M);

for k = 1:K
    C_tot = C_tot + C(:,:,:,k);
end

C_tot_blk = zeros(M*L,M*L);

for m = 1:M
    
    C_tot_blk(1+(m-1)*L:m*L,1+(m-1)*L:m*L) = C_tot(:,:,m);
    
end


%Prepare to save simulation results
SE_level1 = zeros(K,M);
SE_level4 = zeros(K,1);
%% Go through all channel realizations
for i = 1:nbrOfRealizations
    
    
     %Level 4
    
    %Extract channel estimate realizations from all UEs to all APs
    Hhatallj = reshape(H_hat(:,i,:),[M*L K*N]);
    
    if V_d == 1
        
        %Compute MR combining
        V = Hhatallj;
        
    else

        V = ((Hhatallj*F_precoding*F_precoding'*Hhatallj')+C_tot_blk+eyeML)\(Hhatallj*F_precoding);%Dp*
    end
        
 
    %Go through all UEs
    for k = 1:K
        
        v = V(:,(k-1)*N+1:k*N); %Extract combining matrix
        
         %Compute numerator and denominator of instantaneous SINR at Level 4
        numerator = (v'*Hhatallj(:,(k-1)*N+1:k*N)*F_precoding((k-1)*N+1:k*N,(k-1)*N+1:k*N));
        denominator = v'*(Hhatallj*F_precoding*F_precoding'*Hhatallj' + C_tot_blk + eyeML)*v - numerator*numerator';
        SE_level4(k) = SE_level4(k) + prelogFactor*real(log2(det(eye(N) + numerator'/denominator*numerator)))/nbrOfRealizations;
            

    end

       
  %----Level1
    
    %Go through all APs
    for m = 1:M
        
        
        %Extract channel estimate realizations from all UEs to AP m
        Hhatallj1 = reshape(H_Combining(1+(m-1)*L:m*L,i,:),[L K*N]);
        Hallj1 = reshape(H_hat(1+(m-1)*L:m*L,i,:),[L K*N]);

        %Compute MR/L-MMSE combining
        V = Hhatallj1;

        %Go through all UEs
        for k = 1:K
            
            v = V(:,(k-1)*N+1:k*N); %Extract combining matrix
            numerator = v'*Hallj1(:,(k-1)*N+1:k*N);
            denominator = v'*(Hallj1*Hallj1' + C_tot(:,:,m) + eyeL)*v - numerator*numerator';
            
            SE_level1(k,m) = SE_level1(k,m) + prelogFactor*real(log2(det(eye(N)+numerator'/denominator*numerator)))/nbrOfRealizations;
        end
    end
end


%Compute SE for Level 1
SE_l1 = max(SE_level1,[],2);
