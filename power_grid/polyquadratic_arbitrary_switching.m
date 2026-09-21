close all;
clear;
clc;

yalmip('clear');

%% 1. USER PARAMETERS

dataFile = 'C:\Users\Farah\Downloads\ss_matrices_A_conns.mat';
solverName = 'mosek';
numberOfCycles = 40;
strictTol = 1e-8;
feasibilityTol = 1e-7;

load(dataFile,'A_matrices','B');

A_matrices = double(A_matrices);
B = double(B);

N = size(A_matrices,3);
n = size(A_matrices,1);

A = cell(1,N);
for i = 1:N
    A{i} = A_matrices(:,:,i);
end

if n == 10
    x0 = [10;0;-5;0;7;0;-3;0;4;0];
else
    x0 = ones(n,1);
end

fprintf('System: %d modes and %d states.\n',N,n);

%% 2. STABILIZING CONTROLLER FOR EACH MODE

[K,Acl] = buildInitialControllers(A,B);

for i = 1:N
    fprintf('Mode %d: rho(Acl) = %.8f\n', ...
        i,max(abs(eig(Acl{i}))));
end

%% 3. POLY-QUADRATIC STABILITY UNDER ARBITRARY SWITCHING
%
% The following LMI is imposed for every pair i,j:
%
% [ Gi + Gi' - Si,  Gi'*Ai';
%   Ai*Gi,           Sj       ] > 0.
%
% There is one Si and one Gi for every mode i.
% The same Gi is reused for all destination modes j.
%
% Use 1:N to test all modes. To test only a subset, replace 1:N, e.g.
% certificateModes = [1 4 7 3];

certificateModes = 1:N;
AclCertificate = Acl(certificateModes);

[polyFeasible,Slocal,Glocal,epsilonValue] = ...
    solvePolyQuadraticStability( ...
    AclCertificate,n,solverName,strictTol,feasibilityTol);

if ~polyFeasible
    error(['No poly-quadratic certificate was found for the selected ' ...
        'modes. YALMIP/MOSEK message is printed above.']);
end

fprintf('\nPOLY-QUADRATIC STABILITY CERTIFIED.\n');
fprintf('Certified modes: %s\n',mat2str(certificateModes));
fprintf('LMI margin: %.6e\n',epsilonValue);

% P_i = S_i^{-1} are the quadratic Lyapunov matrices used in
% V_i(x) = x''*P_i*x.
Plocal = invertCertificateFamily(Slocal);

% Store the matrices using the original global mode numbers.
P = cell(1,N);
S = cell(1,N);
G = cell(1,N);

for q = 1:numel(certificateModes)
    mode = certificateModes(q);
    P{mode} = Plocal{q};
    S{mode} = Slocal{q};
    G{mode} = Glocal{q};
end

%% 4. NUMERICAL VERIFICATION OF ALL MODE TRANSITIONS
%
% For every i -> j, verify
%       A_i'' P_j A_i - P_i < 0.

[worstEigenvalue,worstPair] = verifyPolyQuadraticDecrease( ...
    Acl,P,certificateModes);

fprintf('\nWorst eigenvalue of A_i'' P_j A_i - P_i: %.6e\n', ...
    worstEigenvalue);
fprintf('Worst transition: mode %d -> mode %d\n', ...
    worstPair(1),worstPair(2));

if worstEigenvalue < 0
    fprintf('All pairwise Lyapunov decrease tests are satisfied.\n');
else
    warning(['The reconstructed P_i matrices have a nonnegative numerical ' ...
        'decrease value. Check solver tolerances and conditioning.']);
end

%% 5. SWITCHING SIMULATION
%
% This sequence is only an example. Since every pair i,j was included in
% the LMIs, any switching sequence using certificateModes is admissible.

requestedBlock = [1 4 7 3];

if all(ismember(requestedBlock,certificateModes))
    switchingBlock = requestedBlock;
else
    switchingBlock = certificateModes;
end

sigma = repmat(switchingBlock,1,numberOfCycles);
x = simulateSwitching(Acl,sigma,x0);

% At x(k), use the Lyapunov function corresponding to the mode that will
% be applied next. Thus, for i = sigma(k) and j = sigma(k+1), the decrease
% comparison is V_j(x(k+1)) < V_i(x(k)).
certificateSchedule = [sigma,sigma(1)];
Vactive = evaluateModeDependentLyapunov(x,P,certificateSchedule);

plotPolyQuadraticSimulation(x,sigma,Vactive);
printLyapunovDecrease(Vactive);

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

        K{i} = -place(A{i},B,poles);
        Acl{i} = A{i} + B*K{i};
    end
end

function [feasible,Svalue,Gvalue,epsilonValue] = ...
    solvePolyQuadraticStability( ...
    Acl,n,solverName,strictTol,feasibilityTol)

    numberOfModes = numel(Acl);

    S = cell(1,numberOfModes);
    G = cell(1,numberOfModes);

    constraints = [];
    traceSum = 0;

    % One positive-definite S_i and one full matrix G_i per mode.
    for i = 1:numberOfModes
        S{i} = sdpvar(n,n,'symmetric');
        G{i} = sdpvar(n,n,'full');

        constraints = [constraints, ...
            S{i} >= strictTol*eye(n)]; %#ok<AGROW>

        traceSum = traceSum + trace(S{i});
    end

    epsilon = sdpvar(1);

    % Normalization removes the homogeneous scaling ambiguity.
    constraints = [constraints, ...
        traceSum == numberOfModes*n, ...
        epsilon >= 0];

    % Exact poly-quadratic LMI:
    %
    % [ G_i + G_i' - S_i,  G_i' A_i';
    %   A_i G_i,            S_j        ] >= epsilon I
    %
    % IMPORTANT: G{i} depends only on i and is reused for every j.
    for i = 1:numberOfModes
        Ai = Acl{i};
        Gi = G{i};

        for j = 1:numberOfModes
            Mij = [ ...
                Gi + Gi' - S{i}, Gi'*Ai'; ...
                Ai*Gi,            S{j}];

            constraints = [constraints, ...
                Mij >= epsilon*eye(2*n)]; %#ok<AGROW>
        end
    end

    options = sdpsettings( ...
        'solver',solverName, ...
        'verbose',1);

    solution = optimize(constraints,-epsilon,options);

    feasible = solution.problem == 0 && ...
        value(epsilon) > feasibilityTol;

    if feasible
        Svalue = cell(1,numberOfModes);
        Gvalue = cell(1,numberOfModes);

        for i = 1:numberOfModes
            Svalue{i} = symmetrize(value(S{i}));
            Gvalue{i} = value(G{i});
        end

        epsilonValue = value(epsilon);
    else
        Svalue = {};
        Gvalue = {};
        epsilonValue = NaN;

        fprintf('\nPoly-quadratic solver message: %s\n', ...
            solution.info);
    end
end

function Plist = invertCertificateFamily(Slist)

    Plist = cell(size(Slist));

    for i = 1:numel(Slist)
        Plist{i} = Slist{i}\eye(size(Slist{i}));
        Plist{i} = symmetrize(Plist{i});
    end
end

function [worstEigenvalue,worstPair] = ...
    verifyPolyQuadraticDecrease(Acl,P,certificateModes)

    worstEigenvalue = -Inf;
    worstPair = [NaN NaN];

    for sourceMode = certificateModes
        Pi = P{sourceMode};
        Ai = Acl{sourceMode};

        for destinationMode = certificateModes
            Pj = P{destinationMode};

            difference = Ai'*Pj*Ai - Pi;
            difference = symmetrize(difference);
            currentEigenvalue = max(real(eig(difference)));

            if currentEigenvalue > worstEigenvalue
                worstEigenvalue = currentEigenvalue;
                worstPair = [sourceMode destinationMode];
            end
        end
    end
end

function x = simulateSwitching(Acl,sigma,x0)

    numberOfSteps = numel(sigma);
    x = zeros(numel(x0),numberOfSteps+1);
    x(:,1) = x0;

    for k = 1:numberOfSteps
        activeMode = sigma(k);
        x(:,k+1) = Acl{activeMode}*x(:,k);
    end
end

function V = evaluateModeDependentLyapunov( ...
    x,P,certificateSchedule)

    numberOfSamples = size(x,2);

    if numel(certificateSchedule) ~= numberOfSamples
        error(['The certificate schedule must contain one mode index ' ...
            'for every state sample.']);
    end

    V = zeros(1,numberOfSamples);

    for k = 1:numberOfSamples
        activeCertificate = certificateSchedule(k);
        Pk = P{activeCertificate};

        if isempty(Pk)
            error('No Lyapunov matrix is available for mode %d.', ...
                activeCertificate);
        end

        V(k) = real(x(:,k)'*Pk*x(:,k));
    end
end

function plotPolyQuadraticSimulation(x,sigma,Vactive)

    numberOfSteps = numel(sigma);
    t = 0:numberOfSteps;

    angleIndices = 1:2:size(x,1);
    frequencyIndices = 2:2:size(x,1);

    figure('Color','w','Name','Poly-quadratic stability');
    layout = tiledlayout(2,2, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    title(layout,'Poly-quadratic Lyapunov simulation');

    nexttile;
    semilogy(t,max(Vactive,1e-14),'LineWidth',1.6);
    grid on;
    xlabel('k');
    ylabel('V_{\sigma(k)}(x(k))');
    title('Active mode-dependent Lyapunov function');

    nexttile;
    plot(t,x(angleIndices,:)','LineWidth',1.1);
    grid on;
    xlabel('k');
    ylabel('\delta_j');
    title('Phase angles');

    nexttile;
    if isempty(frequencyIndices)
        axis off;
    else
        plot(t,x(frequencyIndices,:)','LineWidth',1.1);
        grid on;
        xlabel('k');
        ylabel('\omega_j');
        title('Frequencies');
    end

    nexttile;
    stairs(1:numberOfSteps,sigma,'LineWidth',1.3);
    grid on;
    xlabel('k');
    ylabel('\sigma(k)');
    title('Switching signal');
    xlim([1 numberOfSteps]);
end

function printLyapunovDecrease(V)

    ratio = V(2:end)./max(V(1:end-1),1e-14);
    tolerance = 1e-8;

    fprintf('\nActive poly-quadratic Lyapunov sequence\n');
    fprintf('Initial value: %.8e\n',V(1));
    fprintf('Final value:   %.8e\n',V(end));
    fprintf('Maximum V(k+1)/V(k): %.10f\n',max(ratio));
    fprintf('Observed increases: %d / %d\n', ...
        sum(ratio > 1+tolerance),numel(ratio));
end

function M = symmetrize(M)

    M = (M+M')/2;
end
