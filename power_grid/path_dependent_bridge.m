close all;
clear;
clc;

yalmip('clear');

%% 1. USER PARAMETERS

dataFile = 'C:\Users\Farah\Downloads\ss_matrices_A_conns.mat';
solverName = 'mosek';
numberOfCycles = 40;
strictTol = 1e-8;
feasibilityTol = 1e-8;

% First certificate group used to construct the bridge controller.
modesP1 = [1 4 7];
mode2 = 2;

load(dataFile,'A_matrices','B');

A_matrices = double(A_matrices);
B = double(B);

N = size(A_matrices,3);
n = size(A_matrices,1);
m = size(B,2);

if N < 7
    error('The script requires modes 1, 2, 3, 4, and 7.');
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

[K,Acl] = buildInitialControllers(A,B); %#ok<NASGU>

fprintf('\nIndividual closed-loop matrices:\n');
for i = 1:N
    rhoMode = max(abs(eig(Acl{i})));
    fprintf('Mode %d: rho(Acl) = %.8f',i,rhoMode);

    if rhoMode < 1
        fprintf('  [stable]\n');
    else
        fprintf('  [UNSTABLE]\n');
    end
end

%% 3. COMMON CERTIFICATE FOR MODES 1, 4, AND 7
%
% This certificate is used only to design K2tilde, as in the original
% bridge construction. It is not used as the final path-dependent proof.

[p1Feasible,P1,epsP1,p1Report] = commonCertificateForModes( ...
    Acl,modesP1,n,solverName,strictTol,feasibilityTol);

fprintf('\nCertificate used for the bridge design:\n');
fprintf('Modes: {1,4,7}\n');
printSolverReport(p1Report);

if ~p1Feasible
    error([ ...
        'A common certificate for modes {1,4,7} was not found, so ' ...
        'K2tilde cannot be constructed with P1 fixed.']);
end

fprintf('P1 margin = %.12e\n',epsP1);

%% 4. DESIGN THE BRIDGE CONTROLLER K2TILDE
%
% K2tilde is designed for the physical dynamics of mode 2 while keeping
% P1 fixed. The new closed-loop bridge matrix is
%
%       A2tilde = A2 + B*K2tilde.

[kBridgeFeasible,K2tilde,A2tilde,epsK2tilde,kBridgeReport] = ...
    designK2tilde(A{mode2},B,P1,n,m,solverName,feasibilityTol);

fprintf('\nBridge-controller design:\n');
printSolverReport(kBridgeReport);

if ~kBridgeFeasible
    error('K2tilde could not be constructed with the fixed certificate P1.');
end

fprintf('K2tilde margin   = %.12e\n',epsK2tilde);
fprintf('rho(A2tilde)     = %.12f\n',max(abs(eig(A2tilde))));
fprintf('||K2tilde||_F    = %.12e\n',norm(K2tilde,'fro'));

% Add the bridge mode after the original modes.
mode2tilde = N+1;
AclExtended = Acl;
AclExtended{mode2tilde} = A2tilde;

modeNames = arrayfun(@(i) sprintf('Mode %d',i), ...
    1:N,'UniformOutput',false);
modeNames{mode2tilde} = 'Mode 2 tilde';

%% 5. PERIODIC PATH THAT USES THE K BRIDGE
%
% The active matrix sequence is
%
%   A1, A4, A7, A2tilde, A2, A3, A2, A2tilde, ...
%
% Equivalently, the periodic bridge path is
%
%   1 -> 4 -> 7 -> 2tilde -> 2 -> 3 -> 2 -> 2tilde -> 1.
%
% A separate certificate is associated with every POSITION of this path.
% Therefore, the two occurrences of mode 2 may have different certificates,
% and the two occurrences of mode 2tilde may also have different
% certificates. This is the path-dependent part of the construction.

bridgePathModes = [ ...
    1, ...
    4, ...
    7, ...
    mode2tilde, ...
    2, ...
    3, ...
    2, ...
    mode2tilde];

fprintf('\nBridge path:\n');
printPath(bridgePathModes,modeNames);

%% 6. CHECK THE COMPLETE BRIDGE CYCLE

Acycle = periodicCycleMatrix(AclExtended,bridgePathModes);
rhoCycle = max(abs(eig(Acycle)));
equivalentRate = rhoCycle^(1/numel(bridgePathModes));
cycleGain = norm(Acycle,2);

fprintf('\nComplete bridge-cycle check:\n');
fprintf('rho(Acycle)              = %.12f\n',rhoCycle);
fprintf('Equivalent per-step rate = %.12f\n',equivalentRate);
fprintf('Two-norm cycle gain      = %.12f\n',cycleGain);

if rhoCycle >= 1
    error([ ...
        'The repeated bridge cycle is unstable because rho(Acycle) >= 1. ' ...
        'A strict periodic path-dependent Lyapunov certificate cannot ' ...
        'exist for this exact cycle and these fixed controllers.']);
end

%% 7. PATH-DEPENDENT DILATED LMIs WITH THE BRIDGE
%
% At path position q, the active matrix is A_q and the next certificate is
% S_{q+1}. The imposed condition is
%
% [ G_q + G_q'' - S_q,   G_q'' A_q'';
%   A_q G_q,              S_{q+1}       ] > 0.
%
% Only the eight edges of bridgePathModes are tested. No arbitrary
% pairwise switching condition is imposed.

[bridgeFeasible,Spath,Gpath,epsBridge,pathReport] = ...
    solvePathDependentDilatedLMI( ...
    AclExtended,bridgePathModes,n,solverName,strictTol,feasibilityTol); %#ok<ASGLU>

fprintf('\nPath-dependent bridge-certificate diagnostic:\n');
printSolverReport(pathReport);

if ~bridgeFeasible
    warning([ ...
        'No strict path-dependent certificate satisfying the selected ' ...
        'dilated LMIs was found for the bridge cycle. This does not mean ' ...
        'that no Lyapunov function exists in general. It means only that ' ...
        'this LMI structure, controller design, path, and tolerance did ' ...
        'not produce a certificate.']);
    return;
end

fprintf('\nBRIDGED PATH IS CERTIFIED.\n');
fprintf('Path-dependent LMI margin = %.12e\n',epsBridge);

% The corresponding quadratic functions are
%       V_q(x) = x'' P_q x,   P_q = inv(S_q).
Ppath = invertCertificateFamily(Spath);

%% 8. VERIFY EVERY EDGE DIRECTLY IN THE P VARIABLES

[worstEigenvalue,worstPosition] = verifyPathDecrease( ...
    AclExtended,bridgePathModes,Ppath,modeNames);

fprintf('\nDirect decrease verification:\n');
fprintf('Worst maximum eigenvalue = %.12e\n',worstEigenvalue);
fprintf('Worst path position      = %d\n',worstPosition);

if worstEigenvalue < 0
    fprintf('Every allowed bridge edge satisfies strict Lyapunov decrease.\n');
else
    warning([ ...
        'The recovered inverse matrices do not show strict decrease at ' ...
        'the selected numerical precision. Inspect matrix conditioning ' ...
        'and solver tolerances.']);
end

%% 9. SIMULATE THE BRIDGED PERIODIC PATH

pathLength = numel(bridgePathModes);
positionSequence = repmat(1:pathLength,1,numberOfCycles);
sigmaBridge = bridgePathModes(positionSequence);

xBridge = simulateSwitching(AclExtended,sigmaBridge,x0);

% At x(:,k), use the certificate associated with the current path position.
% After the last matrix of every complete cycle, the certificate returns to
% position 1.
certificateSchedule = [positionSequence,1];
Vactive = evaluatePathDependentLyapunov( ...
    xBridge,Ppath,certificateSchedule);

plotBridgeSimulation( ...
    xBridge,sigmaBridge,Vactive,bridgePathModes,modeNames);

printLyapunovDecrease(Vactive);

%% 10. SAVE THE BRIDGE CONTROLLER AND CERTIFICATES

save('path_dependent_K_bridge_results.mat', ...
    'K2tilde','A2tilde','P1','Spath','Gpath','Ppath', ...
    'bridgePathModes','epsK2tilde','epsBridge','rhoCycle');

fprintf('\nSaved results in path_dependent_K_bridge_results.mat\n');

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
                    poles(2*j-1) = ...
                        radius(j)*exp(1i*angle(j));
                    poles(2*j) = conj(poles(2*j-1));
                end
        end

        % place returns Kplace for A-B*Kplace. We store K=-Kplace because
        % the convention used in this script is Acl=A+B*K.
        K{i} = -place(A{i},B,poles);
        Acl{i} = real(A{i}+B*K{i});
    end
end

function [feasible,Pvalue,epsilonValue,report] = ...
    commonCertificateForModes( ...
    Acl,modes,n,solverName,strictTol,feasibilityTol)

    P = sdpvar(n,n,'symmetric');
    epsilon = sdpvar(1);

    constraints = [ ...
        P >= strictTol*eye(n), ...
        trace(P) == n, ...
        epsilon >= 0];

    for r = 1:numel(modes)
        Ai = Acl{modes(r)};
        AP = Ai*P;

        constraints = [constraints, ...
            [P,AP';AP,P] >= epsilon*eye(2*n)]; %#ok<AGROW>
    end

    options = sdpsettings( ...
        'solver',solverName, ...
        'verbose',1);

    solution = optimize(constraints,-epsilon,options);
    [epsilonComputed,report] = makeSolverReport(solution,epsilon);

    feasible = solution.problem == 0 && ...
        ~isnan(epsilonComputed) && ...
        epsilonComputed > feasibilityTol;

    if feasible
        Pvalue = symmetrize(value(P));
        epsilonValue = epsilonComputed;
    else
        Pvalue = [];
        epsilonValue = NaN;
    end
end

function [feasible,K2tilde,A2tilde,epsilonValue,report] = ...
    designK2tilde( ...
    A2,B,P1,n,m,solverName,feasibilityTol)

    U = sdpvar(m,n,'full');
    epsilon = sdpvar(1);
    APBU = A2*P1+B*U;

    constraints = [ ...
        epsilon >= 0, ...
        [P1,APBU';APBU,P1] >= epsilon*eye(2*n)];

    options = sdpsettings( ...
        'solver',solverName, ...
        'verbose',1);

    solution = optimize(constraints,-epsilon,options);
    [epsilonComputed,report] = makeSolverReport(solution,epsilon);

    feasible = solution.problem == 0 && ...
        ~isnan(epsilonComputed) && ...
        epsilonComputed > feasibilityTol;

    if feasible
        Uvalue = value(U);
        K2tilde = Uvalue/P1;
        A2tilde = real(A2+B*K2tilde);
        epsilonValue = epsilonComputed;
    else
        K2tilde = [];
        A2tilde = [];
        epsilonValue = NaN;
    end
end

function Acycle = periodicCycleMatrix(Acl,pathModes)

    n = size(Acl{pathModes(1)},1);
    Acycle = eye(n);

    for q = 1:numel(pathModes)
        Acycle = Acl{pathModes(q)}*Acycle;
    end
end

function [feasible,Svalue,Gvalue,epsilonValue,report] = ...
    solvePathDependentDilatedLMI( ...
    Acl,pathModes,n,solverName,strictTol,feasibilityTol)

    pathLength = numel(pathModes);

    S = cell(1,pathLength);
    G = cell(1,pathLength);
    constraints = [];
    traceSum = 0;

    for q = 1:pathLength
        S{q} = sdpvar(n,n,'symmetric');
        G{q} = sdpvar(n,n,'full');

        constraints = [constraints, ...
            S{q} >= strictTol*eye(n)]; %#ok<AGROW>

        traceSum = traceSum+trace(S{q});
    end

    epsilon = sdpvar(1);

    constraints = [constraints, ...
        traceSum == pathLength*n, ...
        epsilon >= 0];

    % Only the consecutive edges of the selected periodic bridge path.
    for q = 1:pathLength
        currentMode = pathModes(q);
        nextPosition = mod(q,pathLength)+1;

        Ai = Acl{currentMode};
        Gi = G{q};

        Mq = [ ...
            Gi+Gi'-S{q}, Gi'*Ai'; ...
            Ai*Gi,        S{nextPosition}];

        constraints = [constraints, ...
            Mq >= epsilon*eye(2*n)]; %#ok<AGROW>
    end

    options = sdpsettings( ...
        'solver',solverName, ...
        'verbose',1);

    solution = optimize(constraints,-epsilon,options);
    [epsilonComputed,report] = makeSolverReport(solution,epsilon);

    feasible = solution.problem == 0 && ...
        ~isnan(epsilonComputed) && ...
        epsilonComputed > feasibilityTol;

    if feasible
        Svalue = cell(1,pathLength);
        Gvalue = cell(1,pathLength);

        for q = 1:pathLength
            Svalue{q} = symmetrize(value(S{q}));
            Gvalue{q} = value(G{q});
        end

        epsilonValue = epsilonComputed;
    else
        Svalue = {};
        Gvalue = {};
        epsilonValue = NaN;
    end
end

function Pvalue = invertCertificateFamily(Svalue)

    Pvalue = cell(size(Svalue));

    for q = 1:numel(Svalue)
        Pvalue{q} = Svalue{q}\eye(size(Svalue{q}));
        Pvalue{q} = symmetrize(Pvalue{q});
    end
end

function [worstEigenvalue,worstPosition] = verifyPathDecrease( ...
    Acl,pathModes,Ppath,modeNames)

    pathLength = numel(pathModes);
    worstEigenvalue = -Inf;
    worstPosition = NaN;

    for q = 1:pathLength
        currentMode = pathModes(q);
        nextPosition = mod(q,pathLength)+1;
        nextMode = pathModes(nextPosition);

        Ai = Acl{currentMode};
        Dq = Ai'*Ppath{nextPosition}*Ai-Ppath{q};
        Dq = symmetrize(Dq);
        currentWorst = max(real(eig(Dq)));

        fprintf('q%d: %s -> %s, max eig = %.12e\n', ...
            q,modeNames{currentMode},modeNames{nextMode},currentWorst);

        if currentWorst > worstEigenvalue
            worstEigenvalue = currentWorst;
            worstPosition = q;
        end
    end
end

function x = simulateSwitching(Acl,sigma,x0)

    n = numel(x0);
    numberOfSteps = numel(sigma);

    x = zeros(n,numberOfSteps+1);
    x(:,1) = x0;

    for k = 1:numberOfSteps
        x(:,k+1) = Acl{sigma(k)}*x(:,k);
    end
end

function V = evaluatePathDependentLyapunov( ...
    x,Ppath,certificateSchedule)

    numberOfSamples = size(x,2);

    if numel(certificateSchedule) ~= numberOfSamples
        error('A certificate position is required for every state sample.');
    end

    V = zeros(1,numberOfSamples);

    for k = 1:numberOfSamples
        q = certificateSchedule(k);
        V(k) = real(x(:,k)'*Ppath{q}*x(:,k));
    end
end

function plotBridgeSimulation( ...
    x,sigma,Vactive,pathModes,modeNames)

    numberOfSteps = numel(sigma);
    t = 0:numberOfSteps;

    angleIndices = 1:2:size(x,1);
    frequencyIndices = 2:2:size(x,1);

    figure('Color','w','Name','Path-dependent K bridge');
    layout = tiledlayout(2,2, ...
        'TileSpacing','compact','Padding','compact');
    title(layout,'Path-dependent switching with K2tilde bridge', ...
        'Interpreter','none');

    nexttile;
    semilogy(t,max(Vactive,1e-14),'LineWidth',1.5);
    grid on;
    xlabel('k');
    ylabel('V_{q(k)}(x(k))');
    title('Active path-dependent Lyapunov function');

    nexttile;
    plot(t,x(angleIndices,:)','LineWidth',1.1);
    grid on;
    xlabel('k');
    ylabel('\delta_j');
    title('Phase-angle states');

    nexttile;
    if isempty(frequencyIndices)
        semilogy(t,max(vecnorm(x,2,1),1e-14),'LineWidth',1.5);
        ylabel('||x(k)||_2');
        title('State norm');
    else
        plot(t,x(frequencyIndices,:)','LineWidth',1.1);
        ylabel('\omega_j');
        title('Frequency states');
    end
    grid on;
    xlabel('k');

    nexttile;
    stairs(1:numberOfSteps,sigma,'LineWidth',1.4);
    grid on;
    xlabel('k');
    ylabel('Active mode');
    title('Bridge switching signal');

    activeModes = unique(pathModes,'stable');
    yticks(activeModes);
    yticklabels(modeNames(activeModes));
    xlim([1 numberOfSteps]);
end

function printLyapunovDecrease(V)

    ratio = V(2:end)./max(V(1:end-1),1e-14);
    tolerance = 1e-8;
    numberOfIncreases = sum(ratio > 1+tolerance);

    fprintf('\nSimulation Lyapunov check:\n');
    fprintf('Initial V                    = %.12e\n',V(1));
    fprintf('Final V                      = %.12e\n',V(end));
    fprintf('Maximum V(k+1)/V(k)          = %.12f\n',max(ratio));
    fprintf('Observed numerical increases = %d / %d\n', ...
        numberOfIncreases,numel(ratio));
end

function [epsilonComputed,report] = makeSolverReport(solution,epsilon)

    epsilonComputed = value(epsilon);

    if isempty(epsilonComputed) || ~isfinite(epsilonComputed)
        epsilonComputed = NaN;
    end

    report.problemCode = solution.problem;
    report.explanation = yalmiperror(solution.problem);
    report.solverInfo = solution.info;
    report.epsilon = epsilonComputed;
end

function printSolverReport(report)

    fprintf('Problem code : %d\n',report.problemCode);
    fprintf('Explanation  : %s\n',report.explanation);
    fprintf('Solver info  : %s\n',report.solverInfo);

    if isnan(report.epsilon)
        fprintf('Computed epsilon: unavailable\n');
    else
        fprintf('Computed epsilon: %.12e\n',report.epsilon);
    end
end

function printPath(pathModes,modeNames)

    for q = 1:numel(pathModes)
        fprintf('%s -> ',modeNames{pathModes(q)});
    end

    fprintf('%s\n',modeNames{pathModes(1)});
end

function M = symmetrize(M)
    M = (M+M')/2;
end
