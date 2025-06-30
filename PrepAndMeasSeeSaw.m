function [vstepS,Rho,M]=PrepAndMeasSeeSaw(X,H,Y,B,alpha,beta,dim)
    Rho_SDP = cell(X,1); % sdp variables
    Rho = cell(X,1); % actual states
    Fs = []; % empty list for state constraints
    for i = 1:X
        Rho_SDP{i} = sdpvar(dim,dim,'hermitian','complex');
        Fs = [Fs; Rho_SDP{i}>=0; trace(Rho_SDP{i})==1]; % positivity and unity trace bound
    end
    Sigma_SDP=sdpvar(dim,dim,'hermitian','complex');
    Fs=[Fs;Sigma_SDP>=0.5*(1*Rho_SDP{1}-beta*Rho_SDP{2});]
    Fs=[Fs;Sigma_SDP>=0.5*(1*Rho_SDP{2}-beta*Rho_SDP{1});]
    Fs=[Fs;Sigma_SDP>=0.5*(alpha*Rho_SDP{1}+alpha*Rho_SDP{2});]
    
    Fm=[];
    M_SDP = cell(Y,B);
    M = cell(Y,B);
    for i = 1:Y
        sum = 0;
        R=RandomPOVM(dim,B);
        for j = 1:B
            M{i,j}=R{j};
            M_SDP{i,j} = sdpvar(dim,dim,'hermitian','complex');
            Fm = [Fm; M_SDP{i,j} >= 0;];
            sum = sum + M_SDP{i,j};
        end
        Fm = [Fm; sum == eye(dim);];
    end
    vstepM = 77; vstepS = 777; tol = 0.000000001; % declaring loop parameters 
    while (abs(vstepM-vstepS) >= tol)    
        % sdp for states
        Prob=StratToProb(Rho_SDP,M);
        stepS = real(genOverlap(Prob,alpha,beta)-2*(trace(Sigma_SDP)+2*beta+1)); % the success metric
        % the state optimization step
        diagnostics = optimize([Fs;], -stepS, sdpsettings('solver', 'sdpt3','verbose',0));
        % preparing state for the measurement sdp.
        for x = 1:X
            Rho{x}=value(Rho_SDP{x});
        end
        Sigma = value(Sigma_SDP)
        vstepS = value(stepS)

        % sdp for Alice's measurements
        Prob=StratToProb(Rho,M_SDP);
        stepM = real(genOverlap(Prob,alpha,beta)-2*(trace(Sigma)+2*beta+1)); % the success metric
        % the state optimization step
        diagnostics = optimize([Fm;], -stepM, sdpsettings('solver', 'sdpt3','verbose',0));
        % preparing state for the measurement sdp.
        for y = 1:Y
            for b = 1:B
                M{y,b} = value(M_SDP{y,b});
            end
        end
        % extracting the value of the objective function for the state step
        vstepM = value(stepM)
    end