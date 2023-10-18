%% 最小样本熵的适应度函数
function [ff] = SampleEntropyCost(c,data)
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
dim = 2;   %   dim
tau = 1;   %Downsampling delay time (Users can ignore this if the default value is 1)
for i = 1:K
	x=u(i,:);%
    r = 0.2*std(x);  %   r: Similarity tolerance (usually 0.1*Std(data)~0.25*Std(data))
    fitness(i,:) = SampleEntropy( dim, r, x, tau );
end
[ff] = min(fitness);
end


if nargin < 4, tau = 1; end
if tau > 1, data = downsample(data, tau); end
 
N = length(data);
result = zeros(1,2);
 
for m = dim:dim+1
    Bi = zeros(N-m+1,1);
    dataMat = zeros(N-m+1,m);
    
    % Set the data matrix and construct it as an M-dimensional vector
    for i = 1:N-m+1
        dataMat(i,:) = data(1,i:i+m-1);
    end
    
    % The number of similar modes is calculated using distance
    for j = 1:N-m+1
        % Calculate Chebyshev distance, excluding self-matching cases
        dist = max(abs(dataMat - repmat(dataMat(j,:),N-m+1,1)),[],2);
        % Count the number of dist less than or equal to r
        D = (dist <= r);
        % Self-matching cases are not included
        Bi(j,1) = (sum(D)-1)/(N-m);
    end
 
    % Find the mean of all Bi's
    result(m-dim+1) = sum(Bi)/(N-m+1);
	
end
    % The calculated sample entropy
    sampEn = -log(result(2)/result(1));
	
end



