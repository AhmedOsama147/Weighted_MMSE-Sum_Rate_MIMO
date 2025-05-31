function [U, V, ITER,diff] = myWMMSE(H, tol, SNR, d, sigma2, MAX_ITER)
%% Inputs:
% H : a 4D matrix represents the channel between each user and each BS
% H(a,b,c,d),
% a x b the paths from cth BS to dth user
% c the number of base station
% d the user order at cth BS
% tol : the algorithm thershold 
% SNR : Signal to nosie ration
% d : data vector for each user i at kth BS
% sigma2 : noise power
%% Outputs:
% V : a 4D matrix for transmitter beamforming at each BS to each user where
% U : a 4D matrix for transmitter beamforming at each user from each BS where
% V(a,b,c,d),
% a x b the paths from cth BS to dth user
% c the number of base station
% d the user order at cth BS
%% Function Code:
[N,M,K,IK] = size(H); % Obtain the channel matrix dimensions
% N : number of antennas at the ith user
% M : number of antennas at the kth BS
% K : Number of BSs and number of users
I = IK / K;
P = 10 ^ (SNR/10); %Convert from dB to Linear scale
Pk = P/K; % The equally distributed power among BSs
p_ki = Pk/K; % Equally distributed power among the users at the kth BSs
alpha_ik = ones(K,I); % Equal weights of all users

% Initialize and Generate random V_ik's randoms 
V = 1/sqrt(2) * (randn(M,d,K,I) + 1i * randn(M,d,K,I));

for k = 1:K
    for i = 1:I
        trace_norm = trace(V(:,:,k,i) * V(:,:,k,i)'); % Compute the trace norm (Tr(V V^H))
        V(:,:,k,i) = V(:,:,k,i) * sqrt(p_ki / trace_norm); % Normalize V to satisfy Tr(V V^H)
    end
end

% Initialize the weighted matrix
W = zeros(d,d,K,I);
for k = 1:K
    for i = 1:I
        W(:,:,k,i) = eye(d);
    end
end
ITER = 0; 
while(1)
    ITER = ITER + 1;
    Wold = W;
    U = calculateU(H, V, sigma2); % Function to calculate U
    W = calculateW(H, V, U); % Function to calculate W
    mu_opt = calculateMU(H, U, W, Pk); % Function to calculate mu (Lagrange multiplier)
    %mu_opt = zeros(K,1);
    V = calculateV(H, U, W, alpha_ik, mu_opt, p_ki); % Function to calculate V
    if (WMMSEcondition(W, Wold) <= tol) || (ITER == MAX_ITER) % Stopping condition
        diff = WMMSEcondition(W, Wold);
        break;
    end
end
end