function [coherentx,nonCoherentx] = functionMMSE_interferenceLevels( R_AP,M,K,N,tau_p,Pset)
  
%This function is used to Check the levels of coherent and non-coherent interference levels (used
%for pilot allocation)
%And each AP is equipped with N antennas.
%This is version 1.0 (Last edited: 2022-07-24)


%Prepare to store the results
Ksi = zeros(M,M,K,K);
X_p1 = zeros(M,M,K,K); 
nonCoherent = zeros(K,K);
coherent = zeros(K,K);



%Generate matrix used in this setup

Phi = zeros(N,N,M,K);
Omega = zeros(N,N,M,K);

%Go through all APs          
for m = 1:M
    
    %Go through all UEs
    for k = 1:K
        
        %Compute the UEs indexes that use the same pilot as UE k
        inds = Pset(:,k);
        PsiInv = zeros(N,N);
        
        %Go through all UEs that use the same pilot as UE k 
        for z = 1:length(inds)   
            
            PsiInv = PsiInv +tau_p*R_AP(:,:,m,inds(z)); 
        
        end
            PsiInv = PsiInv + eye(N);
            
            for z = 1:length(inds)
                
                Phi(:,:,m,inds(z)) = PsiInv;
            
            end
            
            Omega(:,:,m,k) = R_AP(:,:,m,k)/PsiInv*R_AP(:,:,m,k);
            
    end
end



% Go through all APs
for m = 1:M
    
    %Go through all UEs
    for k = 1:K
        
        for l=1:K  %Non-coherent interference (i=k')
            
            Ksi(m,m,k,l) = tau_p*trace(R_AP(:,:,m,l)*Omega(:,:,m,k));
            
           
           if any(l==Pset(:,k)) %Coherent interference (If there is pilot contamination)
           
           X_p1(m,m,k,l) = tau_p*trace(R_AP(:,:,m,l)/Phi(:,:,m,k)*R_AP(:,:,m,k));
           

           end
           
        end
        
    end
    
end

%Go through all UEs
for k = 1:K
    
    for l=1:K  %Non-coherent interference (i=k')
    
        nonCoherent(k,l) =trace(Ksi(:,:,k,l)); 
        
        if any(l==Pset(:,k)) 
            
            coherent(k,l)=  abs(trace(X_p1(:,:,k,l)))^2; 
            
        end
        
    end
    
end

coherentx = sum(coherent,2);
nonCoherentx = sum(nonCoherent,2);
              
end
    
