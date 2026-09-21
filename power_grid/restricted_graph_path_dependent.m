close all;
clear;
clc;

yalmip('clear');

%% 1. USER PARAMETERS

dataFile = 'C:\Users\Farah\Downloads\ss_matrices_A_conns.mat';
solverName = 'mosek';
strictTol = 1e-8;
feasibilityTol = 1e-7;
numberOfCycles = 30;

% Two original certified subgroups. Their old Lyapunov certificates are NOT
% used or fixed in this script. All new certificates are free variables.
subgroup1 = [1 4 7];
subgroup2 = [2 3 5 6 8 9];

% The bridge configuration uses the physical dynamics of mode 2 with a new
% feedback gain K2tilde, optimized jointly with the Lyapunov certificates.
bridgeBaseMode = 2;

load(dataFile,'A_matrices','B');

A_matrices = double(A_matrices);
B = double(B);

N = size(A_matrices,3);
n = size(A_matrices,1);
m = size(B,2);

requiredModes = unique([subgroup1 subgroup2 bridgeBaseMode]);
if any(requiredModes > N)
    error('The data file does not contain all requested modes.');
end

if ~isempty(intersect(subgroup1,subgroup2))
    error('The two subgroups must be disjoint.');
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

fprintf('System: %d original modes, %d states, %d inputs.\n',N,n,m);

%% 2. INITIAL MODE CONTROLLERS

[K,Acl] = buildInitialControllers(A,B); %#ok<NASGU>

fprintf('\nIndividual original closed-loop matrices:\n');
for i = 1:N
    rhoMode = max(abs(eig(Acl{i})));
    fprintf('Mode %d: rho(Acl) = %.8f',i,rhoMode);
    if rhoMode < 1
        fprintf('  [stable]\n');
    else
        fprintf('  [UNSTABLE]\n');
    end
end

%% 3. RESTRICTED TRANSITION GRAPH WITH THE BRIDGE
%
% Nodes 1,...,N correspond to the original closed-loop modes.
% Node N+1 corresponds to the new bridge configuration A2+B*K2tilde.
%
% Allowed transitions:
%   - arbitrary transitions inside subgroup 1;
%   - arbitrary transitions inside subgroup 2;
%   - subgroup 1 <-> bridge configuration;
%   - bridge configuration <-> original mode 2.
%
% There are no direct transitions between subgroup 1 and subgroup 2.
% Therefore, every inter-group transition must pass through the bridge.

bridgeNode = N+1;
numberOfNodes = N+1;

allowedEdges = buildRestrictedBridgeEdges( ...
    subgroup1,subgroup2,bridgeNode,bridgeBaseMode);

nodeNames = arrayfun(@(i) sprintf('Mode %d',i), ...
    1:N,'UniformOutput',false);
nodeNames{bridgeNode} = 'Mode 2 tilde';

fprintf('\nRestricted transition graph:\n');
fprintf('Subgroup 1: {%s}\n',strtrim(sprintf('%d ',subgroup1)));
fprintf('Subgroup 2: {%s}\n',strtrim(sprintf('%d ',subgroup2)));
fprintf('Bridge node: %s\n',nodeNames{bridgeNode});
fprintf('Number of allowed directed edges: %d\n',size(allowedEdges,1));

assertNoDirectIntergroupEdges(allowedEdges,subgroup1,subgroup2);
plotRestrictedTransitionGraph( ...
    allowedEdges,numberOfNodes,nodeNames,subgroup1,subgroup2,bridgeNode);

%% 4. JOINT BRIDGE-CONTROLLER AND PATH-DEPENDENT LYAPUNOV DESIGN
%
% All S_i and G_i are optimized from scratch. No old P1, P2, or P22 is
% fixed. The bridge gain is also optimized jointly through
%
%       Ubridge = K2tilde * Gbridge.
%
% For every allowed graph edge i -> j, impose
%
% [ G_i+G_i' - S_i,  (A_i G_i)';
%   A_i G_i,           S_j          ] >= epsilon I.
%
% For the bridge source node, A_i G_i is replaced by
%
%       A_2 G_bridge + B U_bridge,
%
% which keeps the optimization linear. After solving,
%
%       K2tilde = U_bridge / G_bridge.

[feasible,Svalue,Gvalue,K2tilde,A2tilde,epsilonValue,report] = ...
    solveJointRestrictedBridgePDLF( ...
    A,Acl,B,allowedEdges,bridgeNode,bridgeBaseMode, ...
    n,m,solverName,strictTol,feasibilityTol);

fprintf('\nJoint bridge/path-dependent optimization:\n');
printSolverReport(report);

if ~feasible
    error([ ...
        'No feasible restricted-graph path-dependent certificate was ' ...
        'found. This means that the selected controllers, bridge graph, ' ...
        'and LMI structure did not yield a strict certificate; it does ' ...
        'not prove that the switched system is unstable in general.']);
end

fprintf('\nRESTRICTED BRIDGE GRAPH IS CERTIFIED.\n');
fprintf('LMI margin epsilon = %.12e\n',epsilonValue);
fprintf('rho(A2tilde)       = %.12f\n',max(abs(eig(A2tilde))));
fprintf('||K2tilde||_F      = %.12e\n',norm(K2tilde,'fro'));

%% 5. RECOVER AND VERIFY THE LYAPUNOV MATRICES
%
% The graph-dependent quadratic functions are
%
%       V_i(x) = x' P_i x,    P_i = inv(S_i).

Pvalue = invertCertificateFamily(Svalue);

AclExtended = Acl;
AclExtended{bridgeNode} = A2tilde;

[worstEigenvalue,worstEdge,edgeEigenvalues] = verifyAllowedEdges( ...
    AclExtended,allowedEdges,Pvalue,nodeNames);

fprintf('\nDirect verification in P variables:\n');
fprintf('Worst max eigenvalue = %.12e\n',worstEigenvalue);
fprintf('Worst edge           = %s -> %s\n', ...
    nodeNames{worstEdge(1)},nodeNames{worstEdge(2)});

if worstEigenvalue < 0
    fprintf('Every allowed graph transition has strict Lyapunov decrease.\n');
else
    warning([ ...
        'The inverse certificates do not show strict decrease at the ' ...
        'selected numerical precision. Inspect conditioning and solver ' ...
        'tolerances.']);
end

%% 6. SIMULATE ONE ADMISSIBLE BRIDGED SWITCHING SEQUENCE
%
% This is only an illustration. The LMI proof above applies to every
% switching sequence that follows the restricted graph, not only this one.

nodeCycle = [ ...
    1,4,7, ...                 % motion inside subgroup 1
    bridgeNode, ...            % enter the bridge
    2,3,6,8,9,5,2, ...        % motion inside subgroup 2 and return to mode 2
    bridgeNode];               % return through the bridge to subgroup 1

validateNodeCycle(nodeCycle,allowedEdges);

nodeSequence = repmat(nodeCycle,1,numberOfCycles);
x = simulateNodeSequence(AclExtended,nodeSequence,x0);

% At state x(k), use the certificate of the current graph node. After the
% final step, the repeated sequence returns to the first node.
certificateSchedule = [nodeSequence,nodeCycle(1)];
Vactive = evaluateGraphLyapunov(x,Pvalue,certificateSchedule);

plotRestrictedBridgeSimulation( ...
    x,nodeSequence,Vactive,nodeNames,bridgeNode);
printLyapunovDecrease(Vactive);

%% 7. SAVE RESULTS

save('restricted_bridge_path_dependent_results.mat', ...
    'K2tilde','A2tilde','Svalue','Gvalue','Pvalue', ...
    'allowedEdges','subgroup1','subgroup2','bridgeNode', ...
    'bridgeBaseMode','epsilonValue','worstEigenvalue', ...
    'worstEdge','edgeEigenvalues');

fprintf('\nSaved results in restricted_bridge_path_dependent_results.mat\n');

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

        % place returns Kplace for A-B*Kplace. The convention here is
        % Acl=A+B*K, so K=-Kplace.
        K{i} = -place(A{i},B,poles);
        Acl{i} = real(A{i}+B*K{i});
    end
end

function edges = buildRestrictedBridgeEdges( ...
    subgroup1,subgroup2,bridgeNode,bridgeBaseMode)

    edges = [];

    % All directed transitions inside subgroup 1, including self-loops.
    for i = subgroup1
        for j = subgroup1
            edges(end+1,:) = [i j]; %#ok<AGROW>
        end
    end

    % All directed transitions inside subgroup 2, including self-loops.
    for i = subgroup2
        for j = subgroup2
            edges(end+1,:) = [i j]; %#ok<AGROW>
        end
    end

    % Connections between subgroup 1 and the bridge configuration.
    for i = subgroup1
        edges(end+1,:) = [i bridgeNode]; %#ok<AGROW>
        edges(end+1,:) = [bridgeNode i]; %#ok<AGROW>
    end

    % The bridge is connected to subgroup 2 only through original mode 2.
    edges(end+1,:) = [bridgeNode bridgeBaseMode];
    edges(end+1,:) = [bridgeBaseMode bridgeNode];

    edges = unique(edges,'rows','stable');
end

function assertNoDirectIntergroupEdges(edges,subgroup1,subgroup2)

    for e = 1:size(edges,1)
        source = edges(e,1);
        destination = edges(e,2);

        direct12 = ismember(source,subgroup1) && ...
            ismember(destination,subgroup2);
        direct21 = ismember(source,subgroup2) && ...
            ismember(destination,subgroup1);

        if direct12 || direct21
            error('A forbidden direct inter-group transition was created.');
        end
    end
end

function [feasible,Svalue,Gvalue,Kbridge,Abridge,epsilonValue,report] = ...
    solveJointRestrictedBridgePDLF( ...
    A,Acl,B,edges,bridgeNode,bridgeBaseMode, ...
    n,m,solverName,strictTol,feasibilityTol)

    numberOfNodes = bridgeNode;

    S = cell(1,numberOfNodes);
    G = cell(1,numberOfNodes);
    constraints = [];
    traceSum = 0;

    for i = 1:numberOfNodes
        S{i} = sdpvar(n,n,'symmetric');
        G{i} = sdpvar(n,n,'full');

        constraints = [constraints, ...
            S{i} >= strictTol*eye(n)]; %#ok<AGROW>

        traceSum = traceSum+trace(S{i});
    end

    % Ubridge = Kbridge*Gbridge makes the bridge-controller design convex.
    Ubridge = sdpvar(m,n,'full');
    epsilon = sdpvar(1);

    constraints = [constraints, ...
        traceSum == numberOfNodes*n, ...
        epsilon >= 0];

    for e = 1:size(edges,1)
        source = edges(e,1);
        destination = edges(e,2);

        if source == bridgeNode
            AG = A{bridgeBaseMode}*G{bridgeNode}+B*Ubridge;
        else
            AG = Acl{source}*G{source};
        end

        M = [ ...
            G{source}+G{source}'-S{source}, AG'; ...
            AG,                                 S{destination}];

        constraints = [constraints, ...
            M >= epsilon*eye(2*n)]; %#ok<AGROW>
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
        Svalue = cell(1,numberOfNodes);
        Gvalue = cell(1,numberOfNodes);

        for i = 1:numberOfNodes
            Svalue{i} = symmetrize(value(S{i}));
            Gvalue{i} = value(G{i});
        end

        GbridgeValue = Gvalue{bridgeNode};
        UbridgeValue = value(Ubridge);

        if rcond(GbridgeValue) < 1e-12
            warning('Gbridge is poorly conditioned when recovering K2tilde.');
        end

        Kbridge = UbridgeValue/GbridgeValue;
        Abridge = real(A{bridgeBaseMode}+B*Kbridge);
        epsilonValue = epsilonComputed;
    else
        Svalue = {};
        Gvalue = {};
        Kbridge = [];
        Abridge = [];
        epsilonValue = NaN;
    end
end

function Pvalue = invertCertificateFamily(Svalue)

    Pvalue = cell(size(Svalue));

    for i = 1:numel(Svalue)
        Pvalue{i} = Svalue{i}\eye(size(Svalue{i}));
        Pvalue{i} = symmetrize(Pvalue{i});
    end
end

function [worstEigenvalue,worstEdge,edgeEigenvalues] = verifyAllowedEdges( ...
    Acl,edges,Pvalue,nodeNames)

    numberOfEdges = size(edges,1);
    edgeEigenvalues = zeros(numberOfEdges,1);
    worstEigenvalue = -Inf;
    worstEdge = [NaN NaN];

    for e = 1:numberOfEdges
        source = edges(e,1);
        destination = edges(e,2);

        D = Acl{source}'*Pvalue{destination}*Acl{source} ...
            -Pvalue{source};
        D = symmetrize(D);

        edgeEigenvalues(e) = max(real(eig(D)));

        fprintf('%s -> %s: max eig = %.12e\n', ...
            nodeNames{source},nodeNames{destination}, ...
            edgeEigenvalues(e));

        if edgeEigenvalues(e) > worstEigenvalue
            worstEigenvalue = edgeEigenvalues(e);
            worstEdge = [source destination];
        end
    end
end

function validateNodeCycle(nodeCycle,edges)

    nextCycle = [nodeCycle(2:end),nodeCycle(1)];

    for k = 1:numel(nodeCycle)
        edge = [nodeCycle(k),nextCycle(k)];

        if ~ismember(edge,edges,'rows')
            error('The simulation uses a forbidden edge %d -> %d.', ...
                edge(1),edge(2));
        end
    end
end

function x = simulateNodeSequence(Acl,nodeSequence,x0)

    numberOfSteps = numel(nodeSequence);
    x = zeros(numel(x0),numberOfSteps+1);
    x(:,1) = x0;

    for k = 1:numberOfSteps
        x(:,k+1) = Acl{nodeSequence(k)}*x(:,k);
    end
end

function V = evaluateGraphLyapunov(x,Pvalue,certificateSchedule)

    numberOfSamples = size(x,2);

    if numel(certificateSchedule) ~= numberOfSamples
        error('One certificate node is required for every state sample.');
    end

    V = zeros(1,numberOfSamples);

    for k = 1:numberOfSamples
        node = certificateSchedule(k);
        V(k) = real(x(:,k)'*Pvalue{node}*x(:,k));
    end
end

function plotRestrictedTransitionGraph( ...
    edges,numberOfNodes,nodeNames,subgroup1,subgroup2,bridgeNode)

    adjacency = sparse(edges(:,1),edges(:,2),1, ...
        numberOfNodes,numberOfNodes);

    figure('Color','w','Name','Restricted bridge transition graph');
    graphPlot = plot(digraph(adjacency,nodeNames), ...
        'Layout','layered','LineWidth',1.4);

    color1 = [0.0000 0.4470 0.7410];
    color2 = [1.0000 0.0000 0.0000];
    colorBridge = [0.4940 0.1840 0.5560];

    highlight(graphPlot,subgroup1,'NodeColor',color1);
    highlight(graphPlot,subgroup2,'NodeColor',color2);
    highlight(graphPlot,bridgeNode,'NodeColor',colorBridge);

    for e = 1:size(edges,1)
        source = edges(e,1);
        destination = edges(e,2);

        if ismember(source,subgroup1) && ismember(destination,subgroup1)
            edgeColor = color1;
        elseif ismember(source,subgroup2) && ismember(destination,subgroup2)
            edgeColor = color2;
        else
            edgeColor = colorBridge;
        end

        highlight(graphPlot,source,destination,'EdgeColor',edgeColor);
    end

    title('Allowed transitions: inter-group switches only through bridge', ...
        'Interpreter','none');
end

function plotRestrictedBridgeSimulation( ...
    x,nodeSequence,Vactive,nodeNames,bridgeNode)

    numberOfSteps = numel(nodeSequence);
    t = 0:numberOfSteps;
    angleIndices = 1:2:size(x,1);
    frequencyIndices = 2:2:size(x,1);

    figure('Color','w','Name','Restricted bridge PDLF simulation');
    layout = tiledlayout(2,2, ...
        'TileSpacing','compact','Padding','compact');
    title(layout,'Restricted graph with jointly designed K2tilde', ...
        'Interpreter','none');

    nexttile;
    semilogy(t,max(Vactive,1e-14),'LineWidth',1.5);
    grid on;
    xlabel('k');
    ylabel('V_{q(k)}(x(k))');
    title('Active graph-dependent Lyapunov function');

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
    stairs(1:numberOfSteps,nodeSequence,'LineWidth',1.4);
    grid on;
    xlabel('k');
    ylabel('Active graph node');
    title('Admissible switching signal');

    % MATLAB requires tick values to be strictly increasing.
    activeNodes = sort(unique(nodeSequence));
    yticks(activeNodes);

    labels = nodeNames(activeNodes);
    bridgeTick = find(activeNodes == bridgeNode,1);
    if ~isempty(bridgeTick)
        labels{bridgeTick} = 'Mode 2 tilde';
    end
    yticklabels(labels);
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

function M = symmetrize(M)
    M = (M+M')/2;
end
