function [S,Gamma,l,gammaR]=PrepAndMeas(X,H,Y,B,alpha,beta)
    [gammaR,sizeMain,sizeSmall] = RandomMomentMatrix(X+H,Y,B);
    l=length(gammaR)
    gamma = sdpvar(l,l,'hermitian','real');
    con=[];
    con = [con;gamma>=0];
    con=[con;getCons(gammaR,gamma)];
    gammaPrime=gamma(1:sizeMain,1:sizeMain);
    con=[con;gammaPrime>=0]
    LocalGamma=cell(X+H,1);
    Offset=sizeSmall;
    for x = 1:X+H    
        LocalGamma{x} = gamma(sizeMain+Offset*(x-1)+1:sizeMain+Offset*(x),sizeMain+Offset*(x-1)+1:sizeMain+Offset*(x));
    end
    con=[con;LocalGamma{3}>=0.5*(1*LocalGamma{1}-beta*LocalGamma{2});]
    con=[con;LocalGamma{3}>=0.5*(1*LocalGamma{2}-beta*LocalGamma{1});]
    con=[con;LocalGamma{3}>=0.5*(alpha*LocalGamma{1}+alpha*LocalGamma{2});]
    
    Prob = cell(X,Y,B);
    for x = 1:X
        con=[con;gammaPrime(x+1,x+1)==1;LocalGamma{x}>=0;LocalGamma{x}(1,1)==1];
        for y = 1:Y
            Prob{x,y,1}=LocalGamma{x}(1+y);
        end
    end

S=[];
SQ=sdpvar(1,1,'hermitian');
con=[con;LocalGamma{3}(1,1)<=SQ];
%classical bound
S1=real(genOverlap(Prob,alpha,beta)-2*(SQ+2*beta+1));
diagnostics = optimize([con], -S1, sdpsettings('solver', 'mosek'));
Gamma=value(gamma);
S = [S;value(S1)];