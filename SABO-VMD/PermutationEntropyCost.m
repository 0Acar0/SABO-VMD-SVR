%% Fitness function of minimum information entropy
function [ff] = PermutationEntropyCost(c,data)

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
M = 3;  % Embedded dimension
T = 1;  % Delay time
for i = 1:K
    fitness(i,:) = PermutationEntropy(u(i,:),M,T);  % Here the information entropy function is called directly
end
[ff] = min(fitness);
end




%% 排列熵算法
function [pe ,hist] = PermutationEntropy(y,m,t)

%  Calculate the permutation entropy(PE)
%  排列熵算法的提出者：Bandt C，Pompe B. Permutation entropy:a natural complexity measure for time series[J]. Physical Review Letters,2002,88(17):174102.

%  Input:   y: time series;
%           m: order of permuation entropy 嵌入维数
%           t: delay time of permuation entropy,延迟时间

% Output: 
%           pe:    permuation entropy
%           hist:  the histogram for the order distribution
ly = length(y);
permlist = perms(1:m);
[h,~]=size(permlist);
c(1:length(permlist))=0;

 for j=1:ly-t*(m-1)
     [~,iv]=sort(y(j:t:j+t*(m-1)));
     for jj=1:h
         if (abs(permlist(jj,:)-iv))==0
             c(jj) = c(jj) + 1 ;
         end
     end
 end
hist = c;
c=c(c~=0);
p = c/sum(c);
pe = -sum(p .* log(p));
% 归一化
pe=pe/log(factorial(m));
end
