function sum = genOverlap(Prob,alpha,beta)
    sum = (1-alpha) * 0.5 * (Prob{1,1,1} + 1-Prob{2,1,1});
    sum=sum + (1+beta) * 0.5 * Prob{1,2,1} + (alpha+beta) * (1-0.5*(Prob{1,2,1}+Prob{2,2,1}));
    sum=sum + (1+beta) * 0.5 * Prob{2,3,1} + (alpha+beta) * (1-0.5*(Prob{1,3,1}+Prob{2,3,1}));
    sum = real(sum * 2);