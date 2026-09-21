function [feasible,Fopt,Kopt,lambda_opt] = bitsoris_test(A,B,H,w)

n  = size(A,1);
nu = size(B,2);
m  = size(H,1);

K = sdpvar(nu,n,'full');
lambda = sdpvar(1,1);
F = sdpvar(m,m,'full');

constraints = [];
constraints = [constraints, F >= 0];
constraints = [constraints, F*H == H*(A + B*K)];
constraints = [constraints, F*w <= lambda*w];
constraints = [constraints, 1e-6 <= lambda <= 0.999];

% small regularization to avoid undetermined K
objective = lambda + 1e-6*norm(K,1);

ops = sdpsettings('verbose',0,'solver','mosek');
sol = optimize(constraints, objective, ops);

if sol.problem == 0
    Fopt = value(F);
    Kopt = value(K);
    lambda_opt = value(lambda);

    % extra numerical validation
    if all(isfinite(Fopt(:))) && all(isfinite(Kopt(:))) && isfinite(lambda_opt)
        feasible = true;
    else
        feasible = false;
        Fopt = [];
        Kopt = [];
        lambda_opt = [];
    end
else
    feasible = false;
    Fopt = [];
    Kopt = [];
    lambda_opt = [];
    disp(sol.info);
end
end