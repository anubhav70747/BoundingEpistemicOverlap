function [Gamma,sizeMain,sizeSmall] = RandomMomentMatrix(X,Y,B)
    
    [BobMeas,   d] = RandomProjectiveMeasurements(Y, B,4);
    AliceStates = cell(X,1);
    for x = 1:X
           psi = RandomStateVector(d);
           AliceStates{x,1} = psi*psi';
         %   AliceStates{x,1}=RandomDensityMatrix(d);
    end
    dX=RandomPOVM(d,2);
    AliceStates{3,1}=dX{1};
    OperatorList = {eye(d)};
    for x = 1:X
            OperatorList{end+1} = AliceStates{x,1};
    end
    for y = 1:Y
        OperatorList{end+1} = BobMeas{y,1};
    end
    for x = 1:X
        for x1=1:X
            OperatorList{end+1} = AliceStates{x,1}*AliceStates{x1,1};
        end
        for y = 1:Y
            OperatorList{end+1} = AliceStates{x,1}*BobMeas{y,1};
        end       
    end
    sizeMain=length(OperatorList);
%    noOp = length(OperatorList);
%    for i = 2:noOp
%        for j = 2:noOp
%                OperatorList{end+1} = OperatorList{i}*OperatorList{j};
%        end
%    end
    
    for i = 1:length(OperatorList)
        for j = 1:length(OperatorList)
            Gamma(i,j) = trace(OperatorList{i}'*OperatorList{j});
        end
    end

    %for x = 1:7
    %    addGamma=[];
    %    for i = 1:length(OperatorList)
    %        for j = 1:length(OperatorList)
    %        addGamma(i,j) = trace(AliceStates{x,1}*OperatorList{i}'*OperatorList{j});
    %        end
    %    end
    %    Gamma = blkdiag(Gamma,addGamma);
    %end
    LocalOpList={eye(d)};
    for y = 1:Y
            LocalOpList{end+1} = BobMeas{y,1};
    end
    for y=1:Y
        for y1 = 1:Y
            if y~=y1
                LocalOpList{end+1} = BobMeas{y,1}*BobMeas{y1,1};
            end
        end
        for x = 1:X
           % LocalOpList{end+1} = AliceStates{x,1};
          %  LocalOpList{end+1} = AliceStates{x,1}*BobMeas{y,1};
        end
    end

    for x = 1:X
        addGamma=[];
        for i = 1:length(LocalOpList)
            for j = 1:length(LocalOpList)
                addGamma(i,j) = trace(AliceStates{x,1}*LocalOpList{i}'*LocalOpList{j});
            end
        end
        sizeSmall=length(addGamma);
        Gamma = blkdiag(Gamma,addGamma);
    end
