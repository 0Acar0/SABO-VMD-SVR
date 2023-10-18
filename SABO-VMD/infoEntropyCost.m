%% Fitness function of minimum information entropy
function [ff] = infoEntropyCost(c,data)
X = data;
% alpha = 2300;       % moderate bandwidth constraint
alpha = fix(c(1));       % moderate bandwidth constraint
tau = 0;          % noise-tolerance (no strict fidelity enforcement)
K = fix(c(2));              % modes
% K = 10;              % modes
DC = 0;             % no DC part imposed
init = 1;           % initialize omegas uniformly  
tol = 1e-7;     
%--------------- Run actual VMD code---------------------------
[u, u_hat, omega] = VMD(X, alpha, tau, K, DC, init, tol);
for i = 1:K
    fitness(i,:) = entropy(u(i,:));  %Here the information entropy function is called directly
end
[ff] = min(fitness);
end





