close all;
clear;
clc;

yalmip('clear');

%% ============================================================
%  SYSTEM DEFINITION
% =============================================================

A1 = [-1  1;
       1 -1];

B1 = [1;
      0];

A2 = [-1  1;
       1 -1];

B2 = [0;
     -1];

h = 0.05;

% Discretization
A1d = eye(2) + h*A1;
B1d = h*B1;

A2d = eye(2) + h*A2;
B2d = h*B2;

n  = size(A1d,1);
nu = size(B1d,2);


%% ============================================================
%  INITIAL POLYHEDRAL CONSTRAINT SET
%
%       H x <= w
% =============================================================

H0 = [ 1  0;
      -1  0;
       0  1;
       0 -1];

w0 = [1;
      1;
      1;
      1];


%% ============================================================
%  COMMON QUADRATIC LYAPUNOV DESIGN
%
%  Search for:
%
%       Acl_i = Ai + Bi Ki
%
%  sharing the same quadratic Lyapunov function.
% =============================================================

X  = sdpvar(n,n,'symmetric');
Y1 = sdpvar(nu,n,'full');
Y2 = sdpvar(nu,n,'full');

epsi = sdpvar(1,1);

AXBY1 = A1d*X + B1d*Y1;
AXBY2 = A2d*X + B2d*Y2;

L1 = [X,      AXBY1';
      AXBY1,  X];

L2 = [X,      AXBY2';
      AXBY2,  X];

Constraints = [];

Constraints = [Constraints,...
    X >= epsi*eye(n)];

Constraints = [Constraints,...
    L1 >= epsi*eye(2*n)];

Constraints = [Constraints,...
    L2 >= epsi*eye(2*n)];

% Normalization
Constraints = [Constraints,...
    trace(X) == 1];

Constraints = [Constraints,...
    epsi >= 1e-8];


%% Solve

ops = sdpsettings(...
    'verbose',1,...
    'solver','mosek');

sol = optimize(Constraints,-epsi,ops);

if sol.problem ~= 0

    disp(sol.info);
    error('Common quadratic controller design is infeasible.');

end


%% Recover solution

Xv    = value(X);
Y1v   = value(Y1);
Y2v   = value(Y2);
epsiv = value(epsi);

Plyap = inv(Xv);

K1 = Y1v / Xv;
K2 = Y2v / Xv;

Acl1 = A1d + B1d*K1;
Acl2 = A2d + B2d*K2;


%% Display

disp(' ');
disp('=======================================');
disp('COMMON QUADRATIC DESIGN');
disp('=======================================');

disp('epsilon =');
disp(epsiv);

disp('X =');
disp(Xv);

disp('P =');
disp(Plyap);

disp('K1 =');
disp(K1);

disp('K2 =');
disp(K2);

disp('Eigenvalues of Acl1 =');
disp(eig(Acl1));

disp('Eigenvalues of Acl2 =');
disp(eig(Acl2));


%% ============================================================
%  DIRECT LYAPUNOV VERIFICATION
% =============================================================

D1 = Acl1'*Plyap*Acl1 - Plyap;
D2 = Acl2'*Plyap*Acl2 - Plyap;

disp(' ');
disp('Maximum eigenvalue of Lyapunov decrease matrices');

fprintf('Mode 1 : %.6e\n',max(real(eig(D1))));
fprintf('Mode 2 : %.6e\n',max(real(eig(D2))));


%% ============================================================
%  COMMON INVARIANT POLYTOPE
% =============================================================

Pset = Polyhedron('A',H0,'b',w0);

iter    = 0;
maxIter = 100;

while iter < maxIter

    % Check invariance
    inside1 = (Pset >= Acl1*Pset);
    inside2 = (Pset >= Acl2*Pset);

    if inside1 && inside2
        break;
    end

    % Current vertices
    V = Pset.V;

    % Images under the two modes
    Vnext1 = (Acl1*V')';
    Vnext2 = (Acl2*V')';

    % Convex hull of current set and its images
    Vall = [V;
            Vnext1;
            Vnext2];

    Pset = Polyhedron('V',Vall);
    Pset.minHRep();

    iter = iter + 1;

end


fprintf('\nInvariant-set construction stopped after %d iterations.\n',iter);

if iter == maxIter
    warning('Maximum number of iterations reached.');
else
    disp('Common invariant set found.');
end


%% Final H-representation

H = Pset.A;
w = Pset.b;


%% ============================================================
%  PLOT COMMON INVARIANT SET
% =============================================================
%% ============================================================
%  PLOT COMMON QUADRATIC LYAPUNOV ELLIPSOID
% =============================================================

theta = linspace(0,2*pi,400);

% Unit circle
Z = [cos(theta);
     sin(theta)];

% E(P) = {x : x' P x <= 1}
% x = P^(-1/2) z
[V,D] = eig(Plyap);
Pinvhalf = V*diag(1./sqrt(diag(D)))*V';

E = Pinvhalf*Z;

% Images under the two closed-loop modes
E1 = Acl1*E;
E2 = Acl2*E;

figure;
hold on;
grid on;
axis equal;

fill(E(1,:),E(2,:),[0.7 0.7 1],...
    'FaceAlpha',0.25,...
    'EdgeColor','b',...
    'LineWidth',1.5);

fill(E1(1,:),E1(2,:),[1 0.6 0.6],...
    'FaceAlpha',0.25,...
    'EdgeColor','r',...
    'LineWidth',1.5);

fill(E2(1,:),E2(2,:),[0.6 1 0.6],...
    'FaceAlpha',0.25,...
    'EdgeColor',[0 0.5 0],...
    'LineWidth',1.5);

xlabel('x_1');
ylabel('x_2');

title('Common quadratic Lyapunov level set');

legend(...
    '\mathcal{E}(P)',...
    'A_{cl,1}\mathcal{E}(P)',...
    'A_{cl,2}\mathcal{E}(P)',...
    'Location','best',...
    'Interpreter','latex');

set(gca,'FontSize',11);

% Optional:
% exportgraphics(gcf,'quadratic_lyapunov.pdf','ContentType','vector');

figure;
hold on;
grid on;
axis equal;

plot(Pset,...
    'color','blue',...
    'alpha',0.20);

plot(Acl1*Pset,...
    'color','red',...
    'alpha',0.20);

plot(Acl2*Pset,...
    'color','green',...
    'alpha',0.20);

legend(...
    'P',...
    'A_{cl,1}P',...
    'A_{cl,2}P',...
    'Location','best');

title('Common invariant set built with both modes');

xlabel('x_1');
ylabel('x_2');


%% ============================================================
%  NUMERICAL VERTEX CHECK
% =============================================================

Vc = Pset.V';

all_inside1 = true;
all_inside2 = true;

tol = 1e-8;

for k = 1:size(Vc,2)

    x = Vc(:,k);

    if any(H*(Acl1*x) > w + tol)
        all_inside1 = false;
    end

    if any(H*(Acl2*x) > w + tol)
        all_inside2 = false;
    end

end


disp(' ');
disp('=======================================');
disp('NUMERICAL INVARIANCE CHECK');
disp('=======================================');

fprintf('Acl1 P subset P : %d\n',all_inside1);
fprintf('Acl2 P subset P : %d\n',all_inside2);


%% ============================================================
%  BITSORIS TEST - MODE 1
% =============================================================

[feasible1,F1,K1b,lambda1] = ...
    bitsoris_test(A1d,B1d,H,w);


disp(' ');
disp('=======================================');
disp('BITSORIS TEST - MODE 1');
disp('=======================================');

if feasible1

    disp('Bitsoris condition is satisfied.');

    disp('lambda1 =');
    disp(lambda1);

    disp('K1b =');
    disp(K1b);

    disp('F1 =');
    disp(F1);

    Acl1_bits = A1d + B1d*K1b;

    disp('Eigenvalues of new Acl1 =');
    disp(eig(Acl1_bits));

    figure;
    hold on;
    grid on;
    axis equal;

    plot(Pset,...
        'alpha',0.20);

    plot(Acl1_bits*Pset,...
        'alpha',0.20);

    legend(...
        'P',...
        'A_{cl,1}P',...
        'Location','best');

    title(sprintf(...
        'Bitsoris condition - Mode 1, \\lambda = %.4f',...
        lambda1));

    xlabel('x_1');
    ylabel('x_2');

else

    disp('Bitsoris condition is NOT satisfied for mode 1.');

    % Keep the original controller
    Acl1_bits = Acl1;

end


%% ============================================================
%  BITSORIS TEST - MODE 2
% =============================================================

[feasible2,F2,K2b,lambda2] = ...
    bitsoris_test(A2d,B2d,H,w);


disp(' ');
disp('=======================================');
disp('BITSORIS TEST - MODE 2');
disp('=======================================');

if feasible2

    disp('Bitsoris condition is satisfied.');

    disp('lambda2 =');
    disp(lambda2);

    disp('K2b =');
    disp(K2b);

    disp('F2 =');
    disp(F2);

    Acl2_bits = A2d + B2d*K2b;

    disp('Eigenvalues of new Acl2 =');
    disp(eig(Acl2_bits));

    figure;
    hold on;
    grid on;
    axis equal;

    plot(Pset,...
        'alpha',0.20);

    plot(Acl2_bits*Pset,...
        'alpha',0.20);

    legend(...
        'P',...
        'A_{cl,2}P',...
        'Location','best');

    title(sprintf(...
        'Bitsoris condition - Mode 2, \\lambda = %.4f',...
        lambda2));

    xlabel('x_1');
    ylabel('x_2');

else

    disp('Bitsoris condition is NOT satisfied for mode 2.');

    % Keep original controller
    Acl2_bits = Acl2;

end


%% ============================================================
%  VERIFY THE BITSORIS CONTROLLERS NUMERICALLY
% =============================================================

if feasible1

    Vc = Pset.V';

    inside_bits1 = true;

    for k = 1:size(Vc,2)

        x = Vc(:,k);

        if any(H*(Acl1_bits*x) > lambda1*w + 1e-7)
            inside_bits1 = false;
        end

    end

    fprintf(...
        '\nMode 1: Acl1_bits P subset lambda1 P : %d\n',...
        inside_bits1);

end


if feasible2

    Vc = Pset.V';

    inside_bits2 = true;

    for k = 1:size(Vc,2)

        x = Vc(:,k);

        if any(H*(Acl2_bits*x) > lambda2*w + 1e-7)
            inside_bits2 = false;
        end

    end

    fprintf(...
        'Mode 2: Acl2_bits P subset lambda2 P : %d\n',...
        inside_bits2);

end


%% ============================================================
%  LOCAL FUNCTION: BITSORIS TEST
% =============================================================

function [feasible,Fopt,Kopt,lambda_opt] = ...
    bitsoris_test(A,B,H,w)

    n  = size(A,1);
    nu = size(B,2);
    m  = size(H,1);

    % Decision variables
    K      = sdpvar(nu,n,'full');
    F      = sdpvar(m,m,'full');
    lambda = sdpvar(1,1);

    constraints = [];

    % ---------------------------------------------------------
    % IMPORTANT:
    %
    % F >= 0 in the Bitsoris condition means ELEMENTWISE
    % nonnegative.
    %
    % Therefore we use
    %
    %       F(:) >= 0
    %
    % and NOT
    %
    %       F >= 0
    %
    % because YALMIP would interpret the latter as an SDP
    % constraint for a square matrix.
    % ---------------------------------------------------------

    constraints = [constraints,...
        F(:) >= 0];

    % Bitsoris equality:
    %
    %       F H = H(A + BK)
    %
    constraints = [constraints,...
        F*H == H*(A + B*K)];

    % Contractive condition:
    %
    %       F w <= lambda w
    %
    constraints = [constraints,...
        F*w <= lambda*w];

    % Require strict contraction
constraints = [constraints,...
    lambda >= 0];

    % Minimize contraction factor.
    % Small regularization avoids an unnecessarily large K.
    objective = lambda + 1e-6*norm(K,1);

    ops = sdpsettings(...
        'verbose',1,...
        'solver','mosek');

    sol = optimize(constraints,objective,ops);

    if sol.problem == 0

        Fopt       = value(F);
        Kopt       = value(K);
        lambda_opt = value(lambda);

        % Numerical validation
        if all(isfinite(Fopt(:))) && ...
           all(isfinite(Kopt(:))) && ...
           isfinite(lambda_opt)

            feasible = true;

        else

            feasible  = false;
            Fopt      = [];
            Kopt      = [];
            lambda_opt = [];

        end

    else

        feasible  = false;
        Fopt      = [];
        Kopt      = [];
        lambda_opt = [];

        disp('Bitsoris optimization failed:');
        disp(sol.info);

    end

end