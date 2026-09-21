close all;
clear;
clc;

yalmip('clear');

%% 1. USER PARAMETERS

dataFile = 'C:\Users\Farah\Downloads\ss_matrices_A_conns.mat';
solverName = 'mosek';
numberOfCycles = 40;
maxGroupSize = 6;
strictTol = 1e-8;
feasibilityTol = 1e-7;

load(dataFile,'A_matrices','B');

A_matrices = double(A_matrices);
B = double(B);

N = size(A_matrices,3);
n = size(A_matrices,1);
m = size(B,2);

if N < 9
    error('The script requires at least nine modes.');
end

A = cell(1,N);
for i = 1:N
    A{i} = A_matrices(:,:,i);
end

if n == 10
    x0 = [10;0;-5;0;7;0;-3;0;4;0];
else
    x0 = ones(n,1);
end

fprintf('System: %d modes, %d states, %d inputs.\n',N,n,m);

%% 2. INITIAL CONTROLLERS AND CLOSED-LOOP MATRICES

[K,Acl] = buildInitialControllers(A,B);

for i = 1:N
    fprintf('Mode %d: rho(Acl) = %.6f\n', ...
        i,max(abs(eig(Acl{i}))));
end

%% 3. ORIGINAL CERTIFIED GROUPS AND ORIGINAL GRAPH

maxGroupSize = min(maxGroupSize,N);
allCertifiedGroups = findCertifiedGroups( ...
    Acl,N,n,maxGroupSize,solverName,feasibilityTol);

modesP1 = [1 4 7];
modesP2 = [2 3 5 6 8 9];

[P1,epsP1] = exactGroupCertificate(allCertifiedGroups,modesP1);
[P2,epsP2] = exactGroupCertificate(allCertifiedGroups,modesP2);

fprintf('\nOriginal common certificates:\n');
fprintf('P1 for {1,4,7}: epsilon = %.3e\n',epsP1);
fprintf('P2 for {2,3,6,8,9}: epsilon = %.3e\n',epsP2);

nodeNames = arrayfun(@(i) sprintf('Mode %d',i), ...
    1:N,'UniformOutput',false);

Ggroups = plotInitialCertificateGraph( ...
    allCertifiedGroups,N,nodeNames,modesP1,modesP2);

%% 4. DESIGN K2TILDE WITH P1 FIXED

mode2 = 2;
[K2tilde,A2tilde,epsK2tilde] = designK2tilde( ...
    A{mode2},B,P1,n,m,solverName,feasibilityTol);

[P22old,epsP22old] = commonCertificateForMatrices( ...
    {Acl{2},A2tilde},n,solverName,feasibilityTol);

fprintf('\nIntermediate controller:\n');
fprintf('K2tilde with P1 fixed: epsilon = %.3e\n',epsK2tilde);
fprintf('rho(A2tilde) = %.6f\n',max(abs(eig(A2tilde))));
fprintf('Old P22 for {2,2tilde}: epsilon = %.3e\n',epsP22old);

mode2tilde = N + 1;
AclExtended = Acl;
AclExtended{mode2tilde} = A2tilde;
nodeNames{mode2tilde} = 'Mode 2 tilde';

plotAugmentedCertificateGraph( ...
    Ggroups,modesP1,modesP2,mode2,mode2tilde,nodeNames);

%% 5. DILATED PDLF CERTIFICATES FOR DIRECT SWITCHING
%
% Periodic sequence:
%   1 -> 4 -> 7 -> 3 -> 1 -> ...
%
% Certificate states:
%   q1 --A1--> q1 --A4--> q1 --A7--> q1
%   q1 --A3--> q2 --A1--> q1

[directFeasible,Sdirect,epsDirect] = ...
    solveDirectDilatedCertificates( ...
    Acl,n,solverName,strictTol,feasibilityTol);

if directFeasible
    fprintf('\nDIRECT SWITCHING IS CERTIFIED by the dilated LMIs.\n');
    fprintf('Direct LMI margin = %.3e\n',epsDirect);
    Sd1 = Sdirect{1};
    Sd2 = Sdirect{2};
else
    warning(['Direct switching was not certified by the dilated LMIs. ' ...
        'P1 and P2 will be used only for visualization.']);
    Sd1 = P1;
    Sd2 = P2;
end

Qd1 = symmetricInverse(Sd1);
Qd2 = symmetricInverse(Sd2);


[bridgeFeasible,Sbridge,epsBridge] = ...
    solveBridgeDilatedCertificates( ...
    Acl,A2tilde,n,solverName,strictTol,feasibilityTol);

if bridgeFeasible
    fprintf('\nBRIDGED SWITCHING IS CERTIFIED by the dilated LMIs.\n');
    fprintf('Bridge LMI margin = %.3e\n',epsBridge);
    S1 = Sbridge{1};
    S22 = Sbridge{2};
    S2 = Sbridge{3};
else
    warning(['The bridge was not certified by the dilated LMIs. ' ...
        'The old independent certificates will be plotted, but ' ...
        'their curves do not constitute a transition proof.']);
    S1 = P1;
    S22 = P22old;
    S2 = P2;
end

Q1 = symmetricInverse(S1);
Q22 = symmetricInverse(S22);
Q2 = symmetricInverse(S2);

%% 7. DIRECT SWITCHING SIMULATION

directBlock = [1 4 7 3];
sigmaDirect = repmat(directBlock,1,numberOfCycles);
xDirect = simulateSwitching(Acl,sigmaDirect,x0);

% One curve for each certificate, evaluated along the same trajectory.
VdirectAll = evaluateLyapunovFamily(xDirect,{Qd1,Qd2});

% Active certificate sequence for the rigorous switching comparison:
% x(0): q1
% after A1,A4,A7: q1
% after A3: q2
% after the following A1: q1
certDirect = directCertificateSchedule(numberOfCycles);
VdirectActive = evaluateActiveLyapunov( ...
    xDirect,{Qd1,Qd2},certDirect);

plotScenarioWithCertificateCurves( ...
    xDirect,sigmaDirect,nodeNames, ...
    VdirectAll,{'V_1','V_2'},VdirectActive, ...
    'Direct switching: 1 - 4 - 7 - 3');

printLyapunovDecrease(VdirectActive, ...
    'Direct active Lyapunov sequence',directFeasible);

%% 8. BRIDGED SWITCHING SIMULATION

bridgeBlock = [1 4 7 mode2tilde 2 3 2 mode2tilde];
sigmaBridge = repmat(bridgeBlock,1,numberOfCycles);
xBridge = simulateSwitching(AclExtended,sigmaBridge,x0);

% Each certificate separately on the same graph.
VbridgeAll = evaluateLyapunovFamily(xBridge,{Q1,Q22,Q2});

% Active certificate schedule:
% q1, q1, q1, q1, q22, q2, q2, q22, q1, ...
certBridge = bridgeCertificateSchedule(numberOfCycles);
VbridgeActive = evaluateActiveLyapunov( ...
    xBridge,{Q1,Q22,Q2},certBridge);

plotScenarioWithCertificateCurves( ...
    xBridge,sigmaBridge,nodeNames, ...
    VbridgeAll,{'V_1','V_{22}','V_2'},VbridgeActive, ...
    'Bridged switching with K2tilde');

printLyapunovDecrease(VbridgeActive, ...
    'Bridge active Lyapunov sequence',bridgeFeasible);

%% 9. PERIODIC-CYCLE CHECKS

printCycleResult(AclExtended,directBlock,'1 - 4 - 7 - 3');
printCycleResult(AclExtended,bridgeBlock, ...
    '1 - 4 - 7 - 2tilde - 2 - 3 - 2 - 2tilde');

%% ======================== LOCAL FUNCTIONS ==============================

function [K,Acl] = buildInitialControllers(A,B)

    N = numel(A);
    n = size(A{1},1);
    K = cell(1,N);
    Acl = cell(1,N);

    for i = 1:N
        switch mod(i,3)
            case 1
                poles = linspace(0.85,0.97,n);

            case 2
                poles = linspace(0.10,0.35,n);

            otherwise
                if mod(n,2) ~= 0
                    error('Complex pole construction requires even n.');
                end

                radius = linspace(0.75,0.92,n/2);
                angle = linspace(0.25*pi,0.75*pi,n/2);
                poles = zeros(1,n);

                for j = 1:n/2
                    poles(2*j-1) = radius(j)*exp(1i*angle(j));
                    poles(2*j) = conj(poles(2*j-1));
                end
        end

        K{i} = -place(A{i},B,poles);
        Acl{i} = A{i} + B*K{i};
    end
end

function groups = findCertifiedGroups( ...
    Acl,N,n,maxGroupSize,solverName,feasibilityTol)

    groups = struct('modes',{},'P',{},'epsilon',{});
    groupNumber = 0;

    for groupSize = 2:maxGroupSize
        combinations = nchoosek(1:N,groupSize);

        for c = 1:size(combinations,1)
            indices = combinations(c,:);
            matrices = Acl(indices);

            [feasible,Pvalue,epsilonValue] = ...
                tryCommonCertificate( ...
                matrices,n,solverName,feasibilityTol);

            if feasible
                groupNumber = groupNumber + 1;
                groups(groupNumber).modes = indices;
                groups(groupNumber).P = Pvalue;
                groups(groupNumber).epsilon = epsilonValue;
            end
        end
    end

    if isempty(groups)
        error('No certified group was found.');
    end
end

function [Pvalue,epsilonValue] = exactGroupCertificate(groups,targetModes)

    groupID = [];

    for g = 1:numel(groups)
        if isequal(sort(groups(g).modes),sort(targetModes))
            groupID = g;
            break;
        end
    end

    if isempty(groupID)
        error('The exact group {%s} was not found.', ...
            strtrim(sprintf('%d ',targetModes)));
    end

    Pvalue = groups(groupID).P;
    epsilonValue = groups(groupID).epsilon;
end

function [feasible,Pvalue,epsilonValue] = ...
    tryCommonCertificate(matrices,n,solverName,feasibilityTol)

    P = sdpvar(n,n,'symmetric');
    epsilon = sdpvar(1);

    constraints = [ ...
        P >= 1e-8*eye(n), ...
        trace(P) == n, ...
        epsilon >= 0];

    for i = 1:numel(matrices)
        AP = matrices{i}*P;
        constraints = [constraints, ...
            [P,AP';AP,P] >= epsilon*eye(2*n)]; %#ok<AGROW>
    end

    options = sdpsettings('solver',solverName,'verbose',0);
    solution = optimize(constraints,-epsilon,options);

    feasible = solution.problem == 0 && ...
        value(epsilon) > feasibilityTol;

    if feasible
        Pvalue = value(P);
        Pvalue = (Pvalue + Pvalue')/2;
        epsilonValue = value(epsilon);
    else
        Pvalue = [];
        epsilonValue = NaN;
    end
end

function [Pvalue,epsilonValue] = commonCertificateForMatrices( ...
    matrices,n,solverName,feasibilityTol)

    [feasible,Pvalue,epsilonValue] = tryCommonCertificate( ...
        matrices,n,solverName,feasibilityTol);

    if ~feasible
        error('No common certificate was found for the requested matrices.');
    end
end

function [K2tilde,A2tilde,epsilonValue] = designK2tilde( ...
    A2,B,P1,n,m,solverName,feasibilityTol)

    U = sdpvar(m,n,'full');
    epsilon = sdpvar(1);
    APBU = A2*P1 + B*U;

    constraints = [ ...
        epsilon >= 0, ...
        [P1,APBU';APBU,P1] >= epsilon*eye(2*n)];

    options = sdpsettings('solver',solverName,'verbose',0);
    solution = optimize(constraints,-epsilon,options);

    if solution.problem ~= 0 || ...
            value(epsilon) <= feasibilityTol
        error('K2tilde could not be found: %s',solution.info);
    end

    K2tilde = value(U)/P1;
    A2tilde = A2 + B*K2tilde;
    epsilonValue = value(epsilon);
end

function [feasible,Svalue,epsilonValue] = ...
    solveDirectDilatedCertificates( ...
    Acl,n,solverName,strictTol,feasibilityTol)

    S1 = sdpvar(n,n,'symmetric');
    S2 = sdpvar(n,n,'symmetric');
    epsilon = sdpvar(1);

    constraints = [ ...
        S1 >= strictTol*eye(n), ...
        S2 >= strictTol*eye(n), ...
        trace(S1) + trace(S2) == 2*n, ...
        epsilon >= 0];

    % q1 --A1,A4,A7--> q1
    constraints = addDilatedEdge( ...
        constraints,S1,S1,Acl{1},n,epsilon);
    constraints = addDilatedEdge( ...
        constraints,S1,S1,Acl{4},n,epsilon);
    constraints = addDilatedEdge( ...
        constraints,S1,S1,Acl{7},n,epsilon);

    % q1 --A3--> q2
    constraints = addDilatedEdge( ...
        constraints,S1,S2,Acl{3},n,epsilon);

    % q2 --A1--> q1, closing the repeated cycle
    constraints = addDilatedEdge( ...
        constraints,S2,S1,Acl{1},n,epsilon);

    options = sdpsettings('solver',solverName,'verbose',1);
    solution = optimize(constraints,-epsilon,options);

    feasible = solution.problem == 0 && ...
        value(epsilon) > feasibilityTol;

    if feasible
        Svalue = {symmetrize(value(S1)),symmetrize(value(S2))};
        epsilonValue = value(epsilon);
    else
        Svalue = {};
        epsilonValue = NaN;
        fprintf('\nDirect dilated-LMI solver message: %s\n',solution.info);
    end
end

function [feasible,Svalue,epsilonValue] = ...
    solveBridgeDilatedCertificates( ...
    Acl,A2tilde,n,solverName,strictTol,feasibilityTol)

    S1 = sdpvar(n,n,'symmetric');
    S22 = sdpvar(n,n,'symmetric');
    S2 = sdpvar(n,n,'symmetric');
    epsilon = sdpvar(1);

    constraints = [ ...
        S1 >= strictTol*eye(n), ...
        S22 >= strictTol*eye(n), ...
        S2 >= strictTol*eye(n), ...
        trace(S1) + trace(S22) + trace(S2) == 3*n, ...
        epsilon >= 0];

    % q1 --A1,A4,A7--> q1
    constraints = addDilatedEdge( ...
        constraints,S1,S1,Acl{1},n,epsilon);
    constraints = addDilatedEdge( ...
        constraints,S1,S1,Acl{4},n,epsilon);
    constraints = addDilatedEdge( ...
        constraints,S1,S1,Acl{7},n,epsilon);

    % q1 --A2tilde--> q22
    constraints = addDilatedEdge( ...
        constraints,S1,S22,A2tilde,n,epsilon);

    % q22 --A2--> q2
    constraints = addDilatedEdge( ...
        constraints,S22,S2,Acl{2},n,epsilon);

    % q2 --A3--> q2
    constraints = addDilatedEdge( ...
        constraints,S2,S2,Acl{3},n,epsilon);

    % q2 --A2--> q22
    constraints = addDilatedEdge( ...
        constraints,S2,S22,Acl{2},n,epsilon);

    % q22 --A2tilde--> q1
    constraints = addDilatedEdge( ...
        constraints,S22,S1,A2tilde,n,epsilon);

    options = sdpsettings('solver',solverName,'verbose',1);
    solution = optimize(constraints,-epsilon,options);

    feasible = solution.problem == 0 && ...
        value(epsilon) > feasibilityTol;

    if feasible
        Svalue = { ...
            symmetrize(value(S1)), ...
            symmetrize(value(S22)), ...
            symmetrize(value(S2))};
        epsilonValue = value(epsilon);
    else
        Svalue = {};
        epsilonValue = NaN;
        fprintf('\nBridge dilated-LMI solver message: %s\n',solution.info);
    end
end

function constraints = addDilatedEdge( ...
    constraints,Sa,Sb,Aedge,n,epsilon)

    % A separate full slack G is allowed for each graph edge.
    G = sdpvar(n,n,'full');

    M = [ ...
        G + G' - Sa,  G'*Aedge'; ...
        Aedge*G,       Sb];

    constraints = [constraints, ...
        M >= epsilon*eye(2*n)];
end

function Q = symmetricInverse(S)

    Q = S\eye(size(S));
    Q = symmetrize(Q);
end

function M = symmetrize(M)

    M = (M + M')/2;
end

function Ggroups = plotInitialCertificateGraph( ...
    allCertifiedGroups,N,nodeNames,modesP1,modesP2)

    Ggroups = zeros(N,N);

    for g = 1:numel(allCertifiedGroups)
        indices = allCertifiedGroups(g).modes;
        Ggroups(indices,indices) = 1;
    end

    % Keep exactly the original graph type, layout and line width.
    figure('Color','w','Name','Initial certificate graph');
    p = plot(digraph(Ggroups,nodeNames), ...
        'Layout','layered','LineWidth',1.5);

    colorP1 = [0.0000 0.4470 0.7410];   % Blue: certificate P1
    colorP2 = [1.0000 0.0000 0.0000];   % Red: certificate P2

    % Color the nodes certified by P1 and P2.
    highlight(p,modesP1,'NodeColor',colorP1);
    highlight(p,modesP2,'NodeColor',colorP2);

    % Color every directed connection whose two endpoints belong to P1.
    [localSourceP1,localTargetP1] = find(Ggroups(modesP1,modesP1));
    sourceP1 = modesP1(localSourceP1);
    targetP1 = modesP1(localTargetP1);

    if ~isempty(sourceP1)
        highlight(p,sourceP1,targetP1,'EdgeColor',colorP1);
    end

    % Color every directed connection whose two endpoints belong to P2.
    [localSourceP2,localTargetP2] = find(Ggroups(modesP2,modesP2));
    sourceP2 = modesP2(localSourceP2);
    targetP2 = modesP2(localTargetP2);

    if ~isempty(sourceP2)
        highlight(p,sourceP2,targetP2,'EdgeColor',colorP2);
    end

    title('Initial graph of certified mode groups', ...
        'Interpreter','none');
end

function plotAugmentedCertificateGraph( ...
    Ggroups,modesP1,modesP2,mode2,mode2tilde,nodeNames)

    Gnew = zeros(mode2tilde,mode2tilde);
    Gnew(1:mode2tilde-1,1:mode2tilde-1) = Ggroups;

    for i = modesP1
        Gnew(i,mode2tilde) = 1;
        Gnew(mode2tilde,i) = 1;
    end

    Gnew(mode2,mode2tilde) = 1;
    Gnew(mode2tilde,mode2) = 1;
    Gnew = max(Gnew,eye(mode2tilde));

    % Keep exactly the original graph type, layout and line width.
    figure('Color','w','Name','Augmented certificate graph');
    p = plot(digraph(Gnew,nodeNames), ...
        'Layout','layered','LineWidth',1.5);

    colorP1 = [0.0000 0.4470 0.7410];      % Blue: certificate P1
    colorP2 = [1.0000 0.0000 0.0000];      % Red: certificate P2
    colorBridge = [0.4940 0.1840 0.5560];  % Purple: bridge P22

    % Color the nodes. Mode 2 tilde is shown in purple because it is the
    % intermediate bridge mode shared by the two transition certificates.
    highlight(p,modesP1,'NodeColor',colorP1);
    highlight(p,modesP2,'NodeColor',colorP2);
    highlight(p,mode2tilde,'NodeColor',colorBridge);

    % P1: internal connections between modes 1, 4 and 7.
    [localSourceP1,localTargetP1] = find(Gnew(modesP1,modesP1));
    sourceP1 = modesP1(localSourceP1);
    targetP1 = modesP1(localTargetP1);

    if ~isempty(sourceP1)
        highlight(p,sourceP1,targetP1,'EdgeColor',colorP1);
    end

    % P1 also certifies Mode 2 tilde together with the P1 group.
    sourceToBridge = [modesP1, ...
        mode2tilde*ones(1,numel(modesP1))];
    targetToBridge = [mode2tilde*ones(1,numel(modesP1)), ...
        modesP1];
    highlight(p,sourceToBridge,targetToBridge,'EdgeColor',colorP1);

    % P2: all internal directed connections between modes 2, 3, 6, 8 and 9.
    [localSourceP2,localTargetP2] = find(Gnew(modesP2,modesP2));
    sourceP2 = modesP2(localSourceP2);
    targetP2 = modesP2(localTargetP2);

    if ~isempty(sourceP2)
        highlight(p,sourceP2,targetP2,'EdgeColor',colorP2);
    end

    % P22: directed connections between Mode 2 and Mode 2 tilde.
    highlight(p,[mode2 mode2tilde],[mode2tilde mode2], ...
        'EdgeColor',colorBridge);

    title('Graph after adding K2tilde','Interpreter','none');
end

function x = simulateSwitching(Acl,sigma,x0)

    n = numel(x0);
    Nsim = numel(sigma);
    x = zeros(n,Nsim+1);
    x(:,1) = x0;

    for k = 1:Nsim
        activeMode = sigma(k);
        x(:,k+1) = Acl{activeMode}*x(:,k);
    end
end

function Vfamily = evaluateLyapunovFamily(x,Qlist)

    numberOfStates = size(x,2);
    numberOfCertificates = numel(Qlist);
    Vfamily = zeros(numberOfCertificates,numberOfStates);

    for q = 1:numberOfCertificates
        Q = Qlist{q};

        for k = 1:numberOfStates
            Vfamily(q,k) = real(x(:,k)'*Q*x(:,k));
        end
    end
end

function V = evaluateActiveLyapunov(x,Qlist,certificateState)

    numberOfStates = size(x,2);

    if numel(certificateState) ~= numberOfStates
        error('One certificate index is required for each state sample.');
    end

    V = zeros(1,numberOfStates);

    for k = 1:numberOfStates
        q = certificateState(k);
        Q = Qlist{q};
        V(k) = real(x(:,k)'*Q*x(:,k));
    end
end

function schedule = directCertificateSchedule(numberOfCycles)

    % State schedule for repeated mode block [1 4 7 3].
    % Initial state uses q1. Each block ends in q2; the next A1 returns q2->q1.
    modeCertificateAfterStep = repmat([1 1 1 2],1,numberOfCycles);
    schedule = [1,modeCertificateAfterStep];
end

function schedule = bridgeCertificateSchedule(numberOfCycles)

    % Q-list order: {Q1,Q22,Q2}.
    % After [1 4 7 2tilde 2 3 2 2tilde]:
    %       [q1 q1 q1 q22    q2 q2 q22 q1]
    certificateAfterStep = repmat( ...
        [1 1 1 2 3 3 2 1],1,numberOfCycles);
    schedule = [1,certificateAfterStep];
end

function plotScenarioWithCertificateCurves( ...
    x,sigma,nodeNames,Vfamily,Vnames,Vactive,figureTitle)

    Nsim = numel(sigma);
    t = 0:Nsim;
    angleIndices = 1:2:size(x,1);
    frequencyIndices = 2:2:size(x,1);

    figure('Color','w','Name',figureTitle);
    layout = tiledlayout(2,2, ...
        'TileSpacing','compact','Padding','compact');
    title(layout,figureTitle,'Interpreter','none');

    % All individual certificate functions on the same graph.
    nexttile;
    hold on;
    for q = 1:size(Vfamily,1)
        semilogy(t,max(Vfamily(q,:),1e-14),'LineWidth',1.4);
    end
    semilogy(t,max(Vactive,1e-14),'k--','LineWidth',1.8);
    grid on;
    xlabel('k');
    ylabel('V(x(k))');
    title('Lyapunov certificates and active value');
    legend([Vnames,{'Active V'}],'Location','best');

    nexttile;
    plot(t,x(angleIndices,:)','LineWidth',1.2);
    grid on;
    xlabel('k');
    ylabel('\delta_j');
    title('Phase angles');
    legend(makeStateLabels('\delta',numel(angleIndices)), ...
        'Location','best');

    nexttile;
    if isempty(frequencyIndices)
        axis off;
        text(0.5,0.5,'No frequency states', ...
            'HorizontalAlignment','center');
    else
        plot(t,x(frequencyIndices,:)','LineWidth',1.2);
        grid on;
        xlabel('k');
        ylabel('\omega_j');
        title('Frequencies');
        legend(makeStateLabels('\omega',numel(frequencyIndices)), ...
            'Location','best');
    end

    nexttile;
    stairs(1:Nsim,sigma,'LineWidth',1.4);
    grid on;
    xlabel('k');
    ylabel('Active mode');
    title('Switching signal');

    activeModes = sort(unique(sigma));
    yticks(activeModes);
    yticklabels(nodeNames(activeModes));
    xlim([1 Nsim]);
end

function labels = makeStateLabels(symbol,numberOfSignals)

    labels = arrayfun(@(i) sprintf('%s_%d',symbol,i), ...
        1:numberOfSignals,'UniformOutput',false);
end

function printLyapunovDecrease(V,titleText,isCertified)

    ratio = V(2:end)./max(V(1:end-1),1e-14);
    tolerance = 1e-8;
    numberOfIncreases = sum(ratio > 1 + tolerance);

    fprintf('\n%s\n',titleText);
    fprintf('Initial V = %.6e\n',V(1));
    fprintf('Final V   = %.6e\n',V(end));
    fprintf('Maximum V(k+1)/V(k) = %.8f\n',max(ratio));
    fprintf('Observed increases = %d / %d\n', ...
        numberOfIncreases,numel(ratio));

    if isCertified
        fprintf(['The dilated LMIs provide the theoretical ' ...
            'decrease guarantee for the modeled graph edges.\n']);
    else
        fprintf(['No general theoretical decrease guarantee was found; ' ...
            'the plot is only a numerical simulation.\n']);
    end
end

function Acycle = cycleMatrix(Acl,sequence)

    n = size(Acl{sequence(1)},1);
    Acycle = eye(n);

    for k = 1:numel(sequence)
        Acycle = Acl{sequence(k)}*Acycle;
    end
end

function printCycleResult(Acl,sequence,sequenceName)

    Acycle = cycleMatrix(Acl,sequence);
    rho = max(abs(eig(Acycle)));
    gain = norm(Acycle,2);
    equivalentRate = rho^(1/numel(sequence));

    fprintf('\nCycle %s\n',sequenceName);
    fprintf('spectral radius     = %.8f\n',rho);
    fprintf('equivalent per-step = %.8f\n',equivalentRate);
    fprintf('2-norm gain         = %.8f\n',gain);

    if rho < 1
        fprintf('Result: asymptotically stable periodic cycle.\n');
    else
        fprintf('Result: unstable periodic cycle.\n');
    end
end