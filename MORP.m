function result = MORP(ecModel,preResult)
% REQUIREMENT:
% The COBRA toolbox for MATLAB (https://github.com/opencobra/cobratoolbox)
%
% DESCRIPTION:
% This is used for ecModels that are reconstructed by GECKO 3.
% (https://github.com/SysBioChalmers/GECKO)
% The sum of absolute differences between previously and currently 
% simulated enzyme usage fluxes is minimized. While fixing these fluxes
% the total fluxes are subsequently minimized.
%
% INPUTS:
% model: a model that is used in the current simulation
% preResult: structure that contains:
%         * rxns - full list of reactions in the model that was used in the
%                  previous simulation
%         * x - full flux distribution in the previous simulation
%
% OUTPUTS:
% result: structure with:
%         * stat - solution status
%         * obj - objective value (i.e., the sum of absolute differences)
%         * rxns - full list of reactions in the model that is used in the
%                  current simulation
%         * x - full flux distribution in the current simulation
%
% Code partially adapted from linearMOMA.m (COBRA Toolbox 3)

% Match enzyme usage reactions in the previous and current models
[nMets,nRxns] = size(ecModel.S);
commonIdx = ismember(ecModel.rxns,preResult.rxns) & contains(ecModel.rxns,'usage_prot_');
nCommon = sum(commonIdx);
if nCommon == 0
    error('No common enzyme usage reactions in the previous and current models');
else
    tmpMat = speye(nRxns,nRxns);
    commonMat = tmpMat(commonIdx,:); % size = nCommon x nRxns
    commonRxns = ecModel.rxns(commonIdx);
    [~, prePos] = ismember(commonRxns,preResult.rxns);
    commonPreFluxes = preResult.x(prePos);
end


% Solve default FBA problem, which should be feasible.
% If not, then examine the model.
solutionTest = optimizeCbModel(ecModel);

if solutionTest.stat == 1
    % first solving LP to minimize absolute differences
    % Construct the LHS matrix
    % Rows:
    % 1: S * v = 0
    % 2: -v_i + delta_i+ >= -e_i (e represents enzyme usage flux in the previous simulation)
    % 3: v_i  + delta_i- >= e_i
     
    LPproblem.A = [ecModel.S sparse(nMets,2*nCommon);...
        -commonMat speye(nCommon,nCommon) sparse(nCommon,nCommon);...
        commonMat sparse(nCommon,nCommon) speye(nCommon,nCommon)];
    
    % Construct the RHS vector
    LPproblem.b = [ecModel.b; -commonPreFluxes; commonPreFluxes];
    
    % Objective (Minimize delta+ + delta-)
    LPproblem.osense = 1;
    LPproblem.c = [zeros(nRxns,1); ones(2*nCommon,1)];
    
    % Add constraint on delta [0 inf]
    LPproblem.lb = [ecModel.lb; zeros(2*nCommon,1)];
    LPproblem.ub = [ecModel.ub; inf * ones(2*nCommon,1)];
    
    % Construct the constraint direction vector (G for delta, E for others)
    LPproblem.csense = [repmat('E',nMets,1); repmat('G',2*nCommon,1)];     
    
    % Solve
    LPsolution = solveCobraLP(LPproblem);
    
    result_tmp = LPsolution.full(1:nRxns); % flux distribution
    result.obj = LPsolution.obj;
    
    % second solving LP to minimize total fluxes
    ecModel = changeRxnBounds(ecModel,commonRxns,result_tmp(commonIdx),'b');
%     model = changeRxnBounds(model,commonRxns,result_tmp(commonIdx)-abs(result_tmp(commonIdx))*1e-3,'l');
%     model = changeRxnBounds(model,commonRxns,result_tmp(commonIdx)+abs(result_tmp(commonIdx))*1e-3,'u');
    
    LPproblem2.A = [ecModel.S sparse(nMets,2*nRxns);...
        -speye(nRxns,nRxns) speye(nRxns,nRxns) sparse(nRxns,nRxns);...
        speye(nRxns,nRxns) sparse(nRxns,nRxns) speye(nRxns,nRxns)];
    LPproblem2.b = [ecModel.b; zeros(2*nRxns,1)];
    LPproblem2.osense = 1;
    LPproblem2.c = [zeros(nRxns,1); ones(2*nRxns,1)];
    LPproblem2.lb = [ecModel.lb; zeros(2*nRxns,1)];
    LPproblem2.ub = [ecModel.ub; inf * ones(2*nRxns,1)];
    LPproblem2.csense = [repmat('E',nMets,1); repmat('G',2*nRxns,1)];
    LPsolution2 = solveCobraLP(LPproblem2);
    result.x = LPsolution2.full(1:nRxns); % flux distribution
    result.stat = LPsolution2.origStat;
    result.rxns = ecModel.rxns;
else
    error('Default FBA problem of current model is infeasible.');
end

