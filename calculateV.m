function  V = calculateV(H, U, W, alpha_ik, mu_opt, p_ki)
%% Inputs:
% H : a 4D matrix represents the channel between each user and each BS
% H(a,b,c,d),
% a x b the paths from cth BS to dth user
% c the number of base station
% d the user order at cth BS
% alpha_ik : user weight
% mu_opt : Lagrange multiplier for every BS
% U : a 4D matrix for transmitter beamforming at each user from each BS where
%% Outputs:
% V : a 4D matrix for transmitter beamforming at each BS to each user where
% V(a,b,c,d),
% a x b the paths from cth BS to dth user
% c the number of base station
% d the user order at cth BS
%% Function Code:
[N,M,K,IK] = size(H); % Obtain the channel matrix dimensions
% N : number of antennas at the ith user
% M : number of antennas at the kth BS
% K : Number of BSs and number of users
% I : number of users at the kth BS
I = IK / K;
[~,d,~,~] = size(U);
V = zeros(M,d,K,I); 
% d : number of data streams at the ith user at kth BS
% N : number of antennas at the ith user at kth BS
for k = 1:K
    for i = 1:I
        IN = mu_opt(k)*eye(M);
        for j = 1:K
            for l = 1:I
                IN = IN + alpha_ik(j,l)*H(:,:,k,(j-1)*I+l)'*U(:,:,j,l)*W(:,:,j,l)*...
                    U(:,:,j,l)'*H(:,:,k,(j-1)*I+l);
            end
        end
        V(:,:,k,i) = alpha_ik(k,i)* IN \ H(:,:,k,(k-1)*I+i)' * U(:,:,k,i) * W(:,:,k,i);
    end
end
% for k = 1:K
%     for i = 1:I
%         trace_norm = trace(V(:,:,k,i) * V(:,:,k,i)'); % Compute the trace norm (Tr(V V^H))
%         V(:,:,k,i) = V(:,:,k,i) * sqrt(p_ki / trace_norm); % Normalize V to satisfy Tr(V V^H)
%     end
% end
% Check if the matrix contains any NaN values
if any(isnan(V(:)))
    error('V Matrix contains NaN values. Execution paused.');
end
end