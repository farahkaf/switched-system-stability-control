close all;
clear;
clc;

%% Appendix algo of article of blanchini (A separation principle for linear 
%% switching systems and parametrization of all stabilizing controllers)

A1 = [-1  1;
      1 -1];

B1 = [ 1;
      0];

A2 = [-1  1;
      1 -1];

B2 = [ 0;
      -1];

h = 0.05;

% systeme décalé pour imposer un taux de convergence alpha
alpha = 0;
A1dd = A1 + alpha * eye(2);
A2dd = A2 + alpha * eye(2);

% discretisation
A1d = eye(2) + h*A1dd;
B1d = h*B1;

A2d = eye(2) + h*A2dd;
B2d = h*B2;

% Parameters 
lambda  = 0.92;
epsi    = 1e-4;
maxIter = 50;

% Initial symmetric polytope S0
% H = [ 1  0;
%      -1  0;
%       0  1;
%       0 -1];
% 
% w = [1;1;1;1];

H = [ 1  0;
     -1  0;
      0  1;
      0 -1;
      1  1;
     -1 -1;
      1 -1;
     -1  1];

w = [1;1;1;1; 1.5;1.5;1.5;1.5];

S0 = Polyhedron('A',H,'b',w);
Sk = S0;

figure;
hold on;
grid on;
axis equal;
plot(S0,'alpha',0.04,'color','k');
title('iterative construction of S^{(k)}');

for k = 1:maxIter
% current set 
Hk = Sk.A;
hk = Sk.b;

% Normalize second term
Phi = diag(1./hk) * Hk;

disp('Phi_k = ');
disp(Phi);

%% construction of set R

% mode 1
Azu1 = [Phi*A1d, Phi*B1d];
bzu1 = lambda*ones(size(Phi,1),1);

R1 = Polyhedron('A',Azu1,'b',bzu1);

% Mode 2
Azu2 = [Phi*A2d, Phi*B2d];
bzu2 = lambda*ones(size(Phi,1),1);

R2 = Polyhedron('A',Azu2,'b',bzu2);

disp('R1 built successfully');
disp('R2 built successfully');

%% P (projection on x)
P1 = projection(R1,1:2);
P1 = P1.minHRep();

P2 = projection(R2,1:2);
P2 = P2.minHRep();

Pk = intersect(P1,P2);
Pk = Pk.minHRep();

Snext = intersect(Pk,S0);
Snext = Snext.minHRep();

disp('P1 built successfully');
disp('P2 built successfully');
disp('Pk built successfully');
disp('Snext built successfully');

if Snext.isEmptySet()
        error('S^(k+1) is empty. lambda is probably too stringent.');
    end

    plot(Snext,'alpha',0.08);

    fprintf('Iteration %2d: nb vertices = %d\n',k,size(Snext.V,1));

% verification 
   if k == maxIter || isSubsetScaled(Sk,Snext,lambda+epsi)

        disp('Vertices of P1:'); disp(P1.V);
        disp('Vertices of P2:'); disp(P2.V);
        disp('Vertices of Pk:'); disp(Pk.V);

        P1b = intersect(P1,S0); 
        P1b = P1b.minHRep();

        P2b = intersect(P2,S0); 
        P2b = P2b.minHRep();

        Pkb = intersect(Pk,S0);
        Pkb = Pkb.minHRep();

        figure;
        hold on; grid on; axis equal;

        plot(S0,'alpha',0.05,'color','k');
        plot(P1b,'alpha',0.15,'color','b');
        plot(P2b,'alpha',0.15,'color','g');
        plot(Pkb,'alpha',0.20,'color','m');
        plot(Snext,'alpha',0.25,'color','r');

        legend('S0','P1 ∩ S0','P2 ∩ S0','Pk ∩ S0','Snext');
        title(['Final meaningful iteration k = ', num2str(k)]);

        xlim([-1.2 1.2]);
        ylim([-1.2 1.2]);
    end

% test d'arret 
 if isSubsetScaled(Sk,Snext,lambda+epsi)
        fprintf('Converged at iteration %d\n', k);
        Sk = Snext;
        break;
 end

  Sk = Snext;
end
Xset = Sk;

figure;
hold on;
grid on;
axis equal;
plot(S0,'alpha',0.05,'color','b');
plot(Xset,'alpha',0.25,'color','r');
title('Initial set S^(0) and final set X');
legend('S^(0)','X');

%% construction of X, U and P
V = Xset.V;          
X = V';              
nu = size(X,2);

disp('Vertices row-wise V = ');
disp(V);
disp('Matrix X = ');
disp(X);

HX = Xset.A;
hX = Xset.b;

gamma = lambda + epsi;

% mode 1
j = 1;
xj = X(:,j);

disp('vertex x_j = ');
disp(xj);

U1  = zeros(1,nu);
Xi1 = zeros(2,nu);


for j = 1:nu
    xj = X(:,j);

    [u1j, lb, ub, feasible] = admissible_control_in_scaled_set(xj,A1d,B1d,HX,hX,gamma);

    if ~feasible
        error('No admissible control found for mode 1, vertex %d.', j);
    end

    U1(j) = u1j;
    Xi1(:,j) = A1d*xj + B1d*u1j;

    fprintf('vertex %2d: u = %+ .6f,  interval = [%+ .6f , %+ .6f]\n', ...
        j, u1j, lb, ub);
end

disp('U1 = ');
disp(U1);

disp('Xi1 = ');
disp(Xi1);

% mode 2

U2  = zeros(1,nu);
Xi2 = zeros(2,nu);

fprintf('\n===== Step 9: Mode 2 =====\n');

for j = 1:nu
    xj = X(:,j);

    [u2j, lb, ub, feasible] = admissible_control_in_scaled_set(xj, A2d, B2d, HX, hX, gamma);

    if ~feasible
        error('No admissible control found for mode 2, vertex %d.', j);
    end

    U2(j) = u2j;
    Xi2(:,j) = A2d*xj + B2d*u2j;

    fprintf('vertex %2d: u = %+ .6f,  interval = [%+ .6f , %+ .6f]\n', ...
        j, u2j, lb, ub);
end

disp('U2 = ');
disp(U2);

disp('Xi2 = ');
disp(Xi2);

% P1 et P2

P1mat = zeros(nu,nu);
P2mat = zeros(nu,nu);

for j = 1:nu
    % mode 1
    xi = Xi1(:,j);
    pj = decompose_in_X(xi, X, gamma);
    P1mat(:,j) = pj;

    % verification
    err = norm(xi - X*pj);
    fprintf('Mode 1, column %2d: ||xi - X*p|| = %.3e, ||p||_1 = %.6f\n', ...
        j, err, norm(pj,1));
end

for j = 1:nu
    % mode 2
    xi = Xi2(:,j);
    pj = decompose_in_X(xi, X, gamma);
    P2mat(:,j) = pj;

    % verification
    err = norm(xi - X*pj);
    fprintf('Mode 2, column %2d: ||xi - X*p|| = %.3e, ||p||_1 = %.6f\n', ...
        j, err, norm(pj,1));
end

disp('P1 = ');
disp(P1mat);

disp('P2 = ');
disp(P2mat);

fprintf('norm(P1,1) = %.6f\n', norm(P1mat,1));
fprintf('norm(P2,1) = %.6f\n', norm(P2mat,1));

%% construction of Z

% n = size(X,1);
% 
% found = false;
% while ~found
%     Z = randn(nu-n, nu);
%     T = [X; Z];
%     if rank(T) == nu
%         found = true;
%     end
% end
% 
% disp('Z = ');
% disp(Z);
% 
% disp('rank([X;Z]) = ');
% disp(rank([X;Z]));

% construction of Z with QR

n = size(X,1);      % dimension état
nu = size(X,2);     % nombre de sommets

% QR complet
[Q,~] = qr(X');     % Q est nu x nu (orthogonale)

% directions orthogonales (complément)
Z = Q(:,n+1:end)';  % taille (nu-n) x nu

% matrice finale
T = [X; Z];

disp('Z = ');
disp(Z);

disp('size(Z) = ');
disp(size(Z));

disp('rank([X;Z]) = ');
disp(rank(T));
%% construction of V1 and V2
V1mat = Z * P1mat;
V2mat = Z * P2mat;

disp('V1 = ');
disp(V1mat);

disp('V2 = ');
disp(V2mat);

%% controller matrices
M1 = [U1; V1mat] / [X; Z];
M2 = [U2; V2mat] / [X; Z];

% block extraction
K1 = M1(1,1:n);
H1 = M1(1,n+1:end);
G1 = M1(2:end,1:n);
F1 = M1(2:end,n+1:end);

K2 = M2(1,1:n);
H2 = M2(1,n+1:end);
G2 = M2(2:end,1:n);
F2 = M2(2:end,n+1:end);

disp('K1 = '); disp(K1);
disp('H1 = '); disp(H1);
disp('G1 = '); disp(G1);
disp('F1 = '); disp(F1);

disp('K2 = '); disp(K2);
disp('H2 = '); disp(H2);
disp('G2 = '); disp(G2);
disp('F2 = '); disp(F2);
%%  Numerical test: is u(x) linear

%% ----- Mode 1 -----
% K1star = U1 / X;          % best linear fit in least-squares sense
% R1 = U1 - K1star*X;
% 
% fprintf('\n----- Mode 1 -----\n');
% fprintf('K1* = [%.6f   %.6f]\n', K1star(1), K1star(2));
% fprintf('||U1 - K1*X||_2   = %.6e\n', norm(R1,2));
% fprintf('||U1 - K1*X||_inf = %.6e\n', norm(R1,inf));
% fprintf('relative error    = %.6e\n', norm(R1,2)/max(norm(U1,2),1e-12));
% 
% disp('U1 = ');
% disp(U1);
% 
% disp('K1*X = ');
% disp(K1star*X);
% 
% disp('Residual R1 = U1 - K1*X = ');
% disp(R1);
% 
% %% ----- Mode 2 -----
% K2star = U2 / X;
% R2 = U2 - K2star*X;
% 
% fprintf('\n----- Mode 2 -----\n');
% fprintf('K2* = [%.6f   %.6f]\n', K2star(1), K2star(2));
% fprintf('||U2 - K2*X||_2   = %.6e\n', norm(R2,2));
% fprintf('||U2 - K2*X||_inf = %.6e\n', norm(R2,inf));
% fprintf('relative error    = %.6e\n', norm(R2,2)/max(norm(U2,2),1e-12));
% 
% disp('U2 = ');
% disp(U2);
% 
% disp('K2*X = ');
% disp(K2star*X);
% 
% disp('Residual R2 = U2 - K2*X = ');
% disp(R2);
% 
% %%%% Test global u1(x) vs u2(x)
% 
% fprintf('\n===== Testing u1(x) vs u2(x) =====\n');
% 
% Ntest = 50;
% found = false;
% 
% for k = 1:Ntest
% 
%     % point aléatoire dans une boîte
%     x = -1 + 2*rand(2,1);
% 
%     if Xset.contains(x)
% 
%         [u1,~,~,f1] = admissible_control_in_scaled_set(x, A1d, B1d, HX, hX, gamma);
%         [u2,~,~,f2] = admissible_control_in_scaled_set(x, A2d, B2d, HX, hX, gamma);
% 
%         if f1 && f2
%             diff = abs(u1 - u2);
% 
%             fprintf('x = [%.3f %.3f], u1 = %.6f, u2 = %.6f, diff = %.3e\n', ...
%                 x(1), x(2), u1, u2, diff);
% 
%             if diff > 1e-6
%                 fprintf('\n>>> FOUND DIFFERENCE!\n');
%                 found = true;
%                 break;
%             end
%         end
%     end
% end
% 
% if ~found
%     fprintf('\nNo difference found in %d samples (try increasing Ntest).\n', Ntest);
%end

function p = decompose_in_X(xi, X, gamma)

nu = size(X,2);

f = [zeros(2*nu,1)];   % feasibility problem only

Aeq = [X, -X];
beq = xi;

Aineq = [ones(1,2*nu)];
bineq = gamma;

lb = zeros(2*nu,1);
ub = [];

opts = optimoptions('linprog','Display','none');

[z,~,exitflag] = linprog(f,Aineq,bineq,Aeq,beq,lb,ub,opts);

if exitflag <= 0
    error('Unable to decompose xi in X with ||p||_1 <= gamma.');
end

p_plus  = z(1:nu);
p_minus = z(nu+1:end);

p = p_plus - p_minus;
end

%% Closed-loop matrices

Acl1 = [A1d + B1d*K1,   B1d*H1;
        G1,             F1];

Acl2 = [A2d + B2d*K2,   B2d*H2;
        G2,             F2];

disp('Acl1 = ');
disp(Acl1);
disp('Acl2 = ');
disp(Acl2);

disp('eig(Acl1) = ');
disp(eig(Acl1));

disp('eig(Acl2) = ');
disp(eig(Acl2));

%% Simulation closed-loop under switching

Nsim = 80;

x = zeros(2,Nsim+1);
z = zeros(size(F1,1),Nsim+1);
u = zeros(1,Nsim);

x(:,1) = [0.8; -0.6];
z(:,1) = zeros(size(F1,1),1);

sigma = zeros(1,Nsim);

for k = 1:Nsim
    % exemple : commutation alternée
    if mod(k,2)==1
        sigma(k) = 1;
    else
        sigma(k) = 2;
    end

    if sigma(k)==1
        u(k) = H1*z(:,k) + K1*x(:,k);
        x(:,k+1) = A1d*x(:,k) + B1d*u(k);
        z(:,k+1) = F1*z(:,k) + G1*x(:,k);
    else
        u(k) = H2*z(:,k) + K2*x(:,k);
        x(:,k+1) = A2d*x(:,k) + B2d*u(k);
        z(:,k+1) = F2*z(:,k) + G2*x(:,k);
    end
end

figure;
plot(0:Nsim, x(1,:), 'LineWidth', 1.5); hold on;
plot(0:Nsim, x(2,:), 'LineWidth', 1.5);
grid on;
xlabel('k');
ylabel('states');
legend('x_1','x_2');
title('Closed-loop states under switching');

figure;
plot(0:Nsim, z', 'LineWidth', 1.5);
grid on;
xlabel('k');
ylabel('controller state');
title('Controller state z');

figure;
stairs(1:Nsim, u, 'LineWidth', 1.5);
grid on;
xlabel('k');
ylabel('u_k');
title('Control input');

figure;
stairs(1:Nsim, sigma, 'LineWidth', 1.5);
grid on;
xlabel('k');
ylabel('\sigma(k)');
title('Switching signal');