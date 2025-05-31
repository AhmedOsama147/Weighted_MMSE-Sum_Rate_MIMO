function mu_opt = calculateMU(H, U, W, Pk)
%% Inputs:
% H : a 4D matrix represents the channel between each user and each BS
% H(a,b,c,d),
% a x b the paths from cth BS to dth user
% c the number of base station
% d the user order at cth BS
% W : the Weight matrix for each user at every BS
% U : the beamformer at the receiver (RX)
% Pk : Power for each BS
%% Outputs:
% mu_opt : Lagrange multiplier for every BS
%% Function Code:
[N,M,K,IK] = size(H); % Obtain the channel matrix dimensions
% N : number of antennas at the ith user
% M : number of antennas at the kth BS
% K : Number of BSs and number of users
% I : number of users at the kth BS
mu_opt = zeros(K,1); % Initialize each mu value
I = IK/K;
for k = 1:K
    HUW = 0;
    for j = 1:K
        for l = 1:I
            HUW = HUW + H(:,:,k,(k-1)*I+l)'*U(:,:,j,l)*W(:,:,j,l)*U(:,:,j,l)'*H(:,:,k,(k-1)*I+l);
        end
    end
    [D,A] = eig(HUW);
    AI = 0;
    for i = 1:I
        AI = AI + H(:,:,k,(k-1)*I+i)'*U(:,:,k,i)*W(:,:,k,i)*W(:,:,k,i)'*U(:,:,k,i)'*H(:,:,k,(k-1)*I+i);
    end
    PHI = D'*AI*D;
    PHImm = diag(PHI);
    Amm = diag(A);
    mu_low = 0 ; mu_high = 10 ; tol = 1e-4; max_iter = 100;
    mu_opt(k) = bisection(PHImm, Amm, Pk, mu_low, mu_high, tol, max_iter);
end
end