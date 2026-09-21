close all;
clear;
clc;

% example de l'article separation principle
% path dependant lyaponov function 

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
alpha = 2;
A1dd = A1 + alpha * eye(2);
A2dd = A2 + alpha * eye(2);

% discretisation
A1d = eye(2) + h*A1;
B1d = h*B1;

A2d = eye(2) + h*A2;
B2d = h*B2;


H = [ 1  0;
     -1  0;
      0  1;
      0 -1];

w = [1;1;1;1];
P=Polyhedron(H,w);

% lqr 
% Q = eye(2);      
% R = 1;           

Kinit1=-place(A1d,B1d,[0.82;0.8]) % placement de pole 
%Kinit1 = -dlqr(A1d, B1d, Q, R);   % lqr
Acl1=A1d+B1d*Kinit1;
abs(eig(Acl1))

Kinit2=-place(A2d,B2d,[0.85;0.75])  % placement de pole 
%Kinit2 = -dlqr(A2d, B2d, Q, R);    %lqr                                % lqr
Acl2=A2d+B2d*Kinit2;                  
abs(eig(Acl2)) 

figure; plot(P,Acl1*P,Acl2*P)

%% poly-quadratic stability test

n = size(Acl1,1);
epsi = sdpvar(1,1);%1e-6;

S1 = sdpvar(n,n,'symmetric');
S2 = sdpvar(n,n,'symmetric');

G1 = sdpvar(n,n,'full');
G2 = sdpvar(n,n,'full');

Constraints = [];

% positivite
Constraints = [Constraints, S1 >= epsi*eye(n)];
Constraints = [Constraints, S2 >= epsi*eye(n)];

M11 = [G1+G1'-S1,  G1'*Acl1';
       Acl1*G1,    S1];
Constraints = [Constraints, M11 >= epsi*eye(2*n)];

M12 = [G1+G1'-S1,  G1'*Acl1';
       Acl1*G1,    S2];
Constraints = [Constraints, M12 >= epsi*eye(2*n)];

M21 = [G2+G2'-S2,  G2'*Acl2';
       Acl2*G2,    S1];
Constraints = [Constraints, M21 >= epsi*eye(2*n)];

M22 = [G2+G2'-S2,  G2'*Acl2';
       Acl2*G2,    S2];
Constraints = [Constraints, M22 >= epsi*eye(2*n), epsi>=0.01];

ops = sdpsettings('verbose',1,'solver','mosek');
sol = optimize(Constraints,epsi,ops);

if sol.problem == 0
    disp('The closed loop switched system is poly-quadratically stable')

    S1v = value(S1);
    S2v = value(S2);
    G1v = value(G1);
    G2v = value(G2);

    P1 = inv(S1v);
    P2 = inv(S2v);

    disp('S1 ='); disp(S1v)
    disp('S2 ='); disp(S2v)
    disp('G1 ='); disp(G1v)
    disp('G2 ='); disp(G2v)
    disp('P1 ='); disp(P1)
    disp('P2 ='); disp(P2)
else
    disp('theorem is not satisfied')
    disp(sol.info)
end

P = Polyhedron(H,w);
iter = 0;
maxIter = 100;

while ~(P >= Acl1*P) || ~(P >= Acl2*P)
    V = P.V;

    Vnext1 = (Acl1*V')';
    Vnext2 = (Acl2*V')';

    P = Polyhedron([V; Vnext1; Vnext2]);
    P = P.minHRep();

    iter = iter + 1;
    if iter >= maxIter
        break;
    end
end
figure;
hold on;
grid on;
axis equal;
plot(P,'color','blue','alpha',0.2);
plot(Acl1*P,'color','red','alpha',0.2);
plot(Acl2*P,'color','green','alpha',0.2);
legend('P','Acl1 P','Acl2 P');
title('Common invariant set built with both modes');
H = P.A;
w = P.b;


%% Bitsoris condition
[feasible,F1,K1,lambda1] = bitsoris_test(A1d,B1d,H,w);

figure
plot(P,(A1d+B1d*K1)*P)

if feasible
    disp('Bitsoris condition is satisfied.');
    disp('A sufficient condition for positive invariance holds.');
    disp('Matrix F = ');
    disp(F1);
    Acl1=A1d+B1d*K1;
    disp(Acl1);
else
    disp('Bitsoris condition is not satisfied.');
end
[feasible,F2,K2,lambda2] = bitsoris_test(A2d,B2d,H,w);

figure
plot(P,(A2d+B2d*K2)*P)


if feasible
    disp('Bitsoris condition is satisfied.');
    disp('A sufficient condition for positive invariance holds.');
    disp('Matrix F = ');
    disp(F2);
    Acl2=A2d+B2d*K2;
    disp(Acl2);
else
    disp('Bitsoris condition is NOT satisfied.');
end

% numerical check on vertices
Vc = P.V';

all_inside1 = true;
all_inside2 = true;

for i = 1:size(Vc,2)
    x = Vc(:,i);

    if any(H*(Acl1*x) > w + 1e-8)
        all_inside1 = false;
    end

    if any(H*(Acl2*x) > w + 1e-8)
        all_inside2 = false;
    end
end



% Building the set
P = make_set(H,w);
P = make_set(H,w);
% A*P
Pimg1 = image_set(Acl1,P);
Pimg2 = image_set(Acl2,P);

% Plot P and A*P
figure;
hold on;
grid on;
axis equal;
plot_invariance(P,Pimg1);
title('Set P and its image A1 P');
legend('P','A1 P');
figure;
hold on;
grid on;
axis equal;
plot_invariance(P,Pimg2);
title('Set P and its image A2 P');
legend('P','A2 P');


