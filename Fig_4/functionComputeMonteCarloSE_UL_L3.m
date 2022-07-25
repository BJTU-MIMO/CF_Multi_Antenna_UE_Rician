function [SE,signal_MR,Gp_MR_total] = functionComputeMonteCarloSE_UL_L3(H_Combining,HH,HMean,tau_c,tau_p,nbrOfRealizations,N,L,K,M,F_precoding)
%Compute uplink SE for Cell-free mMIMO for the four different receiver
%cooperation levels, using either MR or MMSE/L-MMSE combining
%This is version 1.1 (Last edited: 2022-07-24)


%Compute the prelog factor
prelogFactor = (1-tau_p/tau_c);

%Prepare to store simulation results

SE = zeros(K,1);
signal_MR = zeros(M*N,N,K);
scaling_MR = zeros(M*N,M*N,K);
Gp_MR_total = zeros(M*N,M*N,K,K);
A = zeros(M*N,N,K);
Hw = HMean + HH;
%Go through all channel realizations
for i = 1:nbrOfRealizations
    
    
    %-----------------Levels 1-3
    gp_MR = zeros(M*N,N,K,K);
    
    %Go through all APs
    for m = 1:M
        %Extract channel realizations from all UEs to AP m
        Hallj = reshape(Hw(1+(m-1)*L:m*L,i,:),[L K*N]);
        %Extract channel estimate realizations from all UEs to AP m
        Hhatallj = reshape(H_Combining(1+(m-1)*L:m*L,i,:),[L K*N]);

        
        %Compute MR/L-MMSE combining
        V = Hhatallj;
        
        for k = 1:K
            
            v = V(:,(k-1)*N+1:k*N); %Extract combining matrix
            signal_MR((m-1)*N+1:m*N,:,k) = signal_MR((m-1)*N+1:m*N,:,k) + (v'*Hallj(:,(k-1)*N+1:k*N)*F_precoding((k-1)*N+1:k*N,(k-1)*N+1:k*N))/nbrOfRealizations; 
            scaling_MR((m-1)*N+1:m*N,(m-1)*N+1:m*N,k) = scaling_MR((m-1)*N+1:m*N,(m-1)*N+1:m*N,k) + v'*v/nbrOfRealizations; 
   
            for l = 1:K
                
                gp_MR((m-1)*N+1:m*N,:,k,l) = v'*Hallj(:,(l-1)*N+1:l*N)*F_precoding((k-1)*N+1:k*N,(k-1)*N+1:k*N);
            
            end
        end
    end
    
    
    
    
    for k = 1:K
        for l = 1:K
            

            Gp_MR_total(:,:,k,l) = Gp_MR_total(:,:,k,l) + (gp_MR(:,:,k,l)*gp_MR(:,:,k,l)')/nbrOfRealizations;
              
        end
    end
end


Gp_MR = sum(Gp_MR_total,4);
%------------Calculation of LSFD coefficients
for k = 1:K
    
    b = signal_MR(:,:,k);
    A(:,:,k) = ((Gp_MR(:,:,k)) + scaling_MR(:,:,k))\b;
        
end


%Compute the SE
for k = 1:K
    
    b = signal_MR(:,:,k);
    numerator = A(:,:,k)'*b;
    denominator = A(:,:,k)'*(Gp_MR(:,:,k))*A(:,:,k) - numerator*numerator' + A(:,:,k)'*scaling_MR(:,:,k)*A(:,:,k);
    SE(k) = prelogFactor*real(log2(det(eye(N) + numerator'/denominator*numerator)));
        
end

        
