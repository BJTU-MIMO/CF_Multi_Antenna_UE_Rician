function [C_MMSE_MMSE_Combining] = functionMatrixGeneration(R_Vec,F_Pre,F_precoding,M,K,L,N,tau_p,Pset)

%---This function is used to generate the matrix used in the next
%subsequent calculation 
%This is version 1.1 (Last edited: 2022-07-24)

%Prepare to store the result
Phi = zeros(L*N,L*N,M,K);
Omega = zeros(L*N,L*N,M,K);
C_MMSE = zeros(L*N,L*N,M,K);
C_MMSE_MMSE_Combining = zeros(L,L,M,K); %Matrix for MMSE Combining

 FF = zeros(K*N,K*N);
for k = 1:K
    FF((k-1)*N+1:k*N,(k-1)*N+1:k*N) = F_precoding((k-1)*N+1:k*N,(k-1)*N+1:k*N)*F_precoding((k-1)*N+1:k*N,(k-1)*N+1:k*N)';
end

%Go through all APs          
for m = 1:M
    
    %Go through all UEs
    for k = 1:K
        
        %Compute the UEs indexes that use the same pilot as UE k
        inds = Pset(:,k);
        PsiInv = zeros(L*N,L*N);
        
        %Go through all UEs that use the same pilot as UE k 
        for z = 1:length(inds)   
            
            PsiInv = PsiInv + tau_p*F_Pre(:,:,inds(z))*R_Vec(:,:,m,inds(z))*F_Pre(:,:,inds(z))';%p(inds(z))*

        end
            PsiInv = PsiInv + eye(L*N);

            for z = 1:length(inds)
                
                Phi(:,:,m,inds(z)) = PsiInv;
            
            end
            
            Omega(:,:,m,k) =tau_p*R_Vec(:,:,m,k)*F_Pre(:,:,k)'/PsiInv*F_Pre(:,:,k)*R_Vec(:,:,m,k);% p(k)*
            
    end
end

%Generate estimation error correlation matrices
for k = 1:K
    
    C_MMSE(:,:,:,k) = R_Vec(:,:,:,k) - Omega(:,:,:,k);
   
end

for k = 1:K
    for m = 1:M
        for j = 1:L
            for q = 1:L
                
                for mm = 1:N
                    for p = 1:N
                        
                        C_MMSE_MMSE_Combining(j,q,m,k) = C_MMSE_MMSE_Combining(j,q,m,k) + C_MMSE((p-1)*L+j,(mm-1)*L+q,m,k)*FF((k-1)*N+p,(k-1)*N+mm);
                        
                    end
                end
                
            end
        end
    end
end

                        
                        

            
