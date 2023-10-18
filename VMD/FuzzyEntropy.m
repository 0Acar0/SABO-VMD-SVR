 
%% Fuzzy entropy function
function FuzEn = FuzzyEntropy(data,dim,r,n,tau)
%
% This function calculates fuzzy entropy (FuzEn) of a univariate signal data
%
% Inputs:
%
% data: univariate signal - a vector of size 1 x N (the number of sample points)
% dim: embedding dimension
% r: threshold (it is usually equal to 0.15 of the standard deviation of a signal - because we normalize signals to have a standard deviation of 1, here, r is usually equal to 0.15)
% n: fuzzy power (it is usually equal to 2)
% tau: time lag (it is usually equal to 1)

%
if nargin == 4, tau = 1; end
if nargin == 3, n = 2; tau=1; end
if tau > 1, data = downsample(data, tau); end
 
N = length(data);
result = zeros(1,2);
 
for m = dim:dim+1
    count = zeros(N-m+1,1);
    dataMat = zeros(N-m+1,m);
    
    % Set the data matrix and construct it as an M-dimensional vector
    for i = 1:N-m+1
        dataMat(i,:) = data(1,i:i+m-1);
    end
    % The number of similar modes is calculated using distance
    for j = 1:N-m+1
        % Calculate Chebyshev distance, excluding self-matching cases
        dataMat=dataMat-mean(dataMat,2);
        tempmat=repmat(dataMat(j,:),N-m+1,1);
        dist = max(abs(dataMat - tempmat),[],2);
        D=exp(-(dist.^n)/r);
        count(j) = (sum(D)-1)/(N-m);
    end
    result(m-dim+1) = sum(count)/(N-m+1);
end
    % The calculated fuzzy entropy
    FuzEn = log(result(1)/result(2));
end