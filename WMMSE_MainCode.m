% Weighted MMSE Approach to Distributed Sum-Utility Maximization for a 
% MIMO Interfering Broadcast Channel - 2011
close all; clc; clear;
tol=1e-2; SNR = 5:5:30; 
d = 4;  sigma2 = 1; num = 100;
K = 10; T = 3; R = 2; I = 4; 

MAX_ITER = 100;
SUMITER = 0;
AvreageDiff = 0;

Rate = zeros(K,I);
SumRATE = zeros(num,1);
SR = zeros(length(SNR),1);
for s = 1:length(SNR)
    for n = 1:num 
    H = 1/sqrt(2) * (randn(R,T,K,I*K) + 1i * randn(R,T,K,I*K));
    [U,V, Iter,diff] = myWMMSE(H, tol, SNR(s), d, sigma2, MAX_ITER);
    SUMITER = SUMITER + Iter;
    AvreageDiff = AvreageDiff + diff;
    for k = 1:K
        for i = 1:I
            Int_noise = sigma2*eye(R);
            for j = 1:K
                for l = 1:I
                    if l ~= i && k ~= j
                        Int_noise = Int_noise + H(:,:,j,(k-1)*I+i) * V(:,:,j,l) * ...
                            V(:,:,j,l)' * H(:,:,j,(k-1)*I+i)';
                    end
                end
            end
            Rate(k,i) = log2(det(eye(R) + H(:,:,k,(k-1)*I+i) * V(:,:,k,i) * V(:,:,k,i)' * H(:,:,k,(k-1)*I+i)' / (Int_noise)));
        end
    end
    SumRATE(n) = sum(Rate,"all");
    end
    SR(s) = mean(SumRATE);
end
plot(SNR,real(SR))
grid on