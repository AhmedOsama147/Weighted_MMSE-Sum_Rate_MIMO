function diff = WMMSEcondition(W, Wold)
[~,~,K,I] = size(W);
SumW = 0;
SumWold = 0;
for k = 1:K
    for i = 1:I
        SumW = SumW + log2(det(W(:,:,k,i)));
        SumWold = SumWold + log2(det(Wold(:,:,k,i)));
    end
end
diff = abs(SumW - SumWold);
end