%%Fitness function of minimum envelope entropy
function [ff] = EnvelopeEntropyCost(c,data)
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
	xx= abs(hilbert(u(i,:))); %Minimum envelope entropy calculation formula
	xxx = xx/sum(xx);
    ssum=0;
	for ii = 1:size(xxx,2)
		bb = xxx(1,ii)*log(xxx(1,ii));
        ssum=ssum+bb;
    end
    fitness(i,:) = -ssum;
end
[ff] = min(fitness);
end


