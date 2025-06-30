function [M,d] = RandomProjectiveMeasurements(Y,K,Rank)

d =  K*Rank;

M = cell(Y,K);
for y = 1:Y
    Unitary = RandomUnitary(d);
    for k = 1:K
        M{y,k} = 0;
        for j = 1:Rank
            M{y,k}=M{y,k}+Unitary(:,(k-1)*Rank + j)*Unitary(:,(k-1)*Rank + j)'; 
        end
    end
end