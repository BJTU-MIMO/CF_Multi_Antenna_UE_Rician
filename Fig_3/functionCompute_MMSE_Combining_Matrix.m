function [V_MMSE_Combining] = functionCompute_MMSE_Combining_Matrix(Hhat,C,nbrOfRealizations,L,N,K,M,F_precoding)
%---This function is used to generate the MMSE combining vector used in the next
%subsequent calculation 

%Store identity matrices of different sizes
eyeL = eye(L);
V_MMSE_Combining = zeros(M*L,nbrOfRealizations,K*N);


%Compute sum of all estimation error correlation matrices at every BS
C_tot = sum(C,4);


%Diagonal matrix with transmit powers and its square root
%Dp = diag(pn);

%% Go through all channel realizations
for n = 1:nbrOfRealizations

    %Go through all APs
    for m = 1:M
        
        
        %Extract channel estimate realizations from all UEs to AP l
        Hhatallj = reshape(Hhat(1+(m-1)*L:m*L,n,:),[L K*N]);

        %Compute MMSE combining
        V_MMSE_Combining((m-1)*L+1:m*L,n,:) = ((Hhatallj*F_precoding*F_precoding'*Hhatallj') + C_tot(:,:,m)+eyeL)\(Hhatallj*F_precoding); 
        
    end
end