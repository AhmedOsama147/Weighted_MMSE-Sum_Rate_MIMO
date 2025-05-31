function U = calculateU(H, V, sigma2)
%% Inputs:
% H : a 4D matrix represents the channel between each user and each BS
% H(a,b,c,d),
% a x b the paths from cth BS to dth user
% c the number of base station
% d the user order at cth BS
% V : the beamformer at the transmitter (BS)
% sigma2 : noise power
%% Outputs:
% U : a 4D matrix for transmitter beamforming at each user from each BS where
% U(a,b,c,d),
% a x b the paths from cth BS to dth user
% c the number of base station
% d the user order at cth BS
%% Function Code:
[N,M,K,IK] = size(H); % Obtain the channel matrix dimensions
% N : number of antennas at the ith user
% M : number of antennas at the kth BS
% K : Number of BSs and number of users
% IK : number of all users in the system
I = IK/K;
[~,d,~,~] = size(V);
U = zeros(N,d,K,I); 
% d : number of data streams at the ith user at kth BS
% N : number of antennas at the ith user at kth BS

for k = 1:K % Iterate over all BSs
    for i = 1:I % Iterate over all users in the kth BS
        HV = zeros; % Initialize HV matrix
        for j = 1:K % Iterate over all BSs
            for l = 1:I % Iterate over all users in the lth BS
                HV = HV + H(:,:,j,(k-1)*I+i)*V(:,:,j,l)*V(:,:,j,l)'*H(:,:,j,(k-1)*I+i)';
            end
        end 
        Jki = HV + sigma2*eye(size(HV,1));
        U(:,:,k,i) = Jki\H(:,:,k,(k-1)*I+i)*V(:,:,k,i);
    end
end
% Check if the matrix contains any NaN values
if any(isnan(U(:)))
    error('U Matrix contains NaN values. Execution paused.');
end
end