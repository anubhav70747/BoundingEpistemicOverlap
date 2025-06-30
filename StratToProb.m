function Prob=StratToProb(Rho,M)
X=length(Rho);
[Y,B]=size(M);
Prob=cell(X,Y,B);
for x = 1:X
    for y = 1:Y
        for b = 1:B
            Prob{x,y,b}=trace(Rho{x}*M{y,b});
        end
    end
end