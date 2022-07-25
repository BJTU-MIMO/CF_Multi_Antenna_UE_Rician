function [SE] = functionCompute_Theoretical_SE( R_Vec,F_precoding1,F_precoding_Pilot1,HMean_Withoutphase,M,K,L,N,tau_p,tau_c,Pset)
%---This function is used to computes the theoretical uplink SE for
%MMSE estimator of level3.
%Each AP is equipped with L antennas and each UE is equipped with N antennas.
%This is version 1.1 (Last edited: 2022-07-24)

%Compute the pre-log factor
%assuming only uplink transmission
prelogFactor = (tau_c-tau_p)/(tau_c);

%Prepare to store the result
Phi = zeros(L*N,L*N,M,K);
Omega = zeros(L*N,L*N,M,K);
F_Pre = zeros(L*N,L*N,K);
FF = zeros(N,N,K); 
Gkk = zeros(M*N,N,K);     
Gama_kl = zeros(M*N,M*N,K,K); 
A = zeros(M*N,N,K); % LSFD coefficient matrix
S_k = zeros(M*N,M*N,K); 
E = zeros(N,N,K); 
SE = zeros(K,1); 
F_precoding = zeros(N,N,K);
F_precoding_Pilot = zeros(N,N,K);
X_p1 = zeros(N,N,K,K,M); 
X_p2 = zeros(N,N,K,K,M); 
X_p3_1 = zeros(N,N,K,K,M); 
X_p3_2 = zeros(N,N,K,K,M);
HM = zeros(L,N,M,K);
%%
%----Generate Matrices Applied Below

for k = 1:K
    F_precoding(:,:,k) = F_precoding1((k-1)*N+1:k*N,(k-1)*N+1:k*N);
    F_precoding_Pilot(:,:,k) = F_precoding_Pilot1((k-1)*N+1:k*N,(k-1)*N+1:k*N);
    FF(:,:,k) = F_precoding(:,:,k)*F_precoding(:,:,k)'; 
    F_Pre(:,:,k) = kron(F_precoding_Pilot(:,:,k).',eye(L));  
    
end

for m = 1:M
    
    for k = 1:K
        
        %Compute the UEs indexes that use the same pilot as UE k
        inds = Pset(:,k);
        PsiInv = zeros(L*N,L*N);
        
        %Go through all UEs that use the same pilot as UE k
        for z = 1:length(inds)
            
            PsiInv = PsiInv + tau_p*F_Pre(:,:,inds(z))*R_Vec(:,:,m,inds(z))*F_Pre(:,:,inds(z))';
            
        end
        
        PsiInv = PsiInv + eye(L*N);
        
        for z = 1:length(inds)
            
            Phi(:,:,m,inds(z)) = PsiInv;
            
        end
        
        Omega(:,:,m,k) = tau_p*R_Vec(:,:,m,k)*F_Pre(:,:,k)'/PsiInv*F_Pre(:,:,k)*R_Vec(:,:,m,k);
        
    end
end
%%

for m = 1:M
    for k = 1:K
        HM(:,:,m,k) = HMean_Withoutphase((m-1)*L+1:m*L,1,(k-1)*N+1:k*N);
    end
end


for k = 1:K
    
    for m = 1:M
        
        for n = 1:N
            for nn = 1:N
                
                Gkk((m-1)*N+n,nn,k) = trace(Omega((nn-1)*L+1:nn*L,(n-1)*L+1:n*L,m,k) + HM(:,nn,m,k)*HM(:,n,m,k)'); 
                
                for l = 1:K
                    
                    for i = 1:N
                        for ii = 1:N  
                            
                            X_p1(n,nn,k,l,m) = X_p1(n,nn,k,l,m) + FF(ii,i,l)*trace((R_Vec((ii-1)*L+1:ii*L,(i-1)*L+1:i*L,m,l) + HM(:,ii,m,l)*HM(:,i,m,l)')*(Omega((nn-1)*L+1:nn*L,(n-1)*L+1:n*L,m,k) + HM(:,nn,m,k)*HM(:,n,m,k)'));
                            
                        end
                    end
                    
                    
                    if any(l == Pset(:,k))
                        phi_sqrt_k = Omega(:,:,m,k)^(1/2);
                        phi_sqrt_l = Omega(:,:,m,l)^(1/2);
                        for i = 1:N
                            for ii = 1:N
                                X_p2_2 = 0;
                                X_p2_3 = 0;
                                for p1 = 1:N
                                    
                                    if any(l == k)
                                        X_p2_2 = X_p2_2 + FF(ii,i,l)*(HM(:,n,m,k)'*HM(:,ii,m,l)*trace(phi_sqrt_l((p1-1)*L+1:p1*L,(i-1)*L+1:i*L)*phi_sqrt_k((nn-1)*L+1:nn*L,(p1-1)*L+1:p1*L))...
                                        +trace(phi_sqrt_k((p1-1)*L+1:p1*L,(n-1)*L+1:n*L)*phi_sqrt_l((ii-1)*L+1:ii*L,(p1-1)*L+1:p1*L))*HM(:,i,m,l)'*HM(:,nn,m,k)...
                                        + HM(:,n,m,k)'*phi_sqrt_l((ii-1)*L+1:ii*L,(p1-1)*L+1:p1*L)*phi_sqrt_l((p1-1)*L+1:p1*L,(i-1)*L+1:i*L)*HM(:,nn,m,k)  ...
                                        +trace(phi_sqrt_k((p1-1)*L+1:p1*L,(n-1)*L+1:n*L)*HM(:,ii,m,l)*HM(:,i,m,l)'*phi_sqrt_k((nn-1)*L+1:nn*L,(p1-1)*L+1:p1*L)));
                                    else
                                        X_p2_2 = X_p2_2 + FF(ii,i,l)*( HM(:,n,m,k)'*phi_sqrt_l((ii-1)*L+1:ii*L,(p1-1)*L+1:p1*L)*phi_sqrt_l((p1-1)*L+1:p1*L,(i-1)*L+1:i*L)*HM(:,nn,m,k)  ...
                                        +trace(phi_sqrt_k((p1-1)*L+1:p1*L,(n-1)*L+1:n*L)*HM(:,ii,m,l)*HM(:,i,m,l)'*phi_sqrt_k((nn-1)*L+1:nn*L,(p1-1)*L+1:p1*L)));
                                    end
                                    for p2 = 1:N
                                         X_p2_3 = X_p2_3 + FF(ii,i,l)*( trace(phi_sqrt_k((p1-1)*L+1:p1*L,(n-1)*L+1:n*L)*phi_sqrt_l((ii-1)*L+1:ii*L,(p2-1)*L+1:p2*L)*phi_sqrt_l((p2-1)*L+1:p2*L,(i-1)*L+1:i*L)*phi_sqrt_k((nn-1)*L+1:nn*L,(p1-1)*L+1:p1*L))...
                                             + trace(phi_sqrt_k((p1-1)*L+1:p1*L,(n-1)*L+1:n*L)*phi_sqrt_l((ii-1)*L+1:ii*L,(p1-1)*L+1:p1*L))*trace(phi_sqrt_k((nn-1)*L+1:nn*L,(p2-1)*L+1:p2*L)*phi_sqrt_l((p2-1)*L+1:p2*L,(i-1)*L+1:i*L)));
                                    end
                                end
                                 
                                    X_p2(n,nn,k,l,m) = X_p2(n,nn,k,l,m) +  X_p2_2 +X_p2_3...
                                     + FF(ii,i,l)*( HM(:,n,m,k)'*HM(:,ii,m,l)*HM(:,i,m,l)'*HM(:,nn,m,k)+trace((Omega((nn-1)*L+1:nn*L,(n-1)*L+1:n*L,m,k)...
                                     +HM(:,nn,m,k)*HM(:,n,m,k)')*(R_Vec((ii-1)*L+1:ii*L,(i-1)*L+1:i*L,m,l)-Omega((ii-1)*L+1:ii*L,(i-1)*L+1:i*L,m,l)))  );
                               
                            end
                        end
                    for ll = 1:N
                        X_p3_1(n,nn,k,l,m) = X_p3_1(n,nn,k,l,m) + trace(phi_sqrt_k((ll-1)*L+1:ll*L,(n-1)*L+1:n*L)*phi_sqrt_l((nn-1)*L+1:nn*L,(ll-1)*L+1:ll*L));
                    end
                    if any(l == k)
                        X_p3_1(n,nn,k,l,m) = X_p3_1(n,nn,k,l,m) + HM(:,n,m,k)'*HM(:,nn,m,l);
                    
                    end
                    
                    end
                end
            end
        end
    end
end

for k = 1:K
    for m = 1:M
        for mm = 1:M
            for l = 1:K
                
                if any(l == Pset(:,k))
                    if any(mm == m)
                        
                        Gama_kl((m-1)*N+1:m*N,(m-1)*N+1:m*N,k,l) = X_p2(:,:,k,l,m);
                    else
                        
                        Gama_kl((m-1)*N+1:m*N,(mm-1)*N+1:mm*N,k,l) = X_p3_1(:,:,k,l,m)*FF(:,:,l)*X_p3_1(:,:,k,l,mm)';
                        
                    end
                else
                    
                    if any(mm == m)
                        
                        Gama_kl((m-1)*N+1:m*N,(m-1)*N+1:m*N,k,l) = X_p1(:,:,k,l,m);
                        
                    end
                    
                end
            end
            S_k((m-1)*N+1:m*N,(m-1)*N+1:m*N,k) = Gkk((m-1)*N+1:m*N,:,k);
        end
    end
end

Gama_kl_total = sum(Gama_kl,4);
for k = 1:K
    
    A(:,:,k) = ((Gama_kl_total(:,:,k)) + S_k(:,:,k))\(Gkk(:,:,k)*F_precoding(:,:,k));
    E(:,:,k) = eye(N) - (Gkk(:,:,k)*F_precoding(:,:,k))'*A(:,:,k);
    
end




%%

%Compute the SE
for k = 1:K
    
    numerator = A(:,:,k)'*Gkk(:,:,k)*F_precoding(:,:,k);
    denominator = A(:,:,k)'*(Gama_kl_total(:,:,k))*A(:,:,k) - numerator*numerator' + A(:,:,k)'*S_k(:,:,k)*A(:,:,k);
    SE(k) = prelogFactor*real(log2(det(eye(N) + numerator'/denominator*numerator)));
    
end

