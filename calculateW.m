function W = calculateW(H, V, U)
%% Inputs:
% H : a 4D matrix represents the channel between each user and each BS
% H(a,b,c,d),
% a x b the paths from cth BS to dth user
% c the number of base station
% d the user order at cth BS
% V : the beamformer at the transmitter (BS)
% U : the beamformer at the receiver (RX)
%% Outputs:
% W : a 4D weight matrix for each user at each BS where
% W(a,a,c,d),
% a x a the matrix dimensions
% c the number of base station
% d the user order at cth BS
%% Function Code:
[N,M,K,IK] = size(H); % Obtain the channel matrix dimensions
% N : number of antennas at the ith user
% M : number of antennas at the kth BS
% K : Number of BSs and number of users
% I : number of users at the kth BS
I = IK/K;
[~,d,~,~] = size(V); % Obtain the data streams
W = zeros(d,d,K,I);
for k = 1:K
    for i = 1:I
        invW = eye(d) - U(:,:,k,i)'*H(:,:,k,(k-1)*I+i)*V(:,:,k,i);
        W(:,:,k,i) = inv(invW);
    end
end
% Check if the matrix contains any NaN values
if any(isnan(W(:)))
    error('W Matrix contains NaN values. Execution paused.');
end
end