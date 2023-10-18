%% sample entropy function
function sampEn = SampleEntropy( dim, r, data, tau )
Note: This sample entropy function is a modification of Kijoon Lee's work
The author of % sample entropy algorithm: Richman J s, Moorman J R. Physiological time-seriesanalysis using approximate entropy and sample entropy[J. American Journal of  Physiology Heart &. Circula-tory Physiology, 2000,278 (6):2039-2049.
% Calculates the sample entropy of the given time series data
% sample entropy is conceptually similar to approximate entropy, but with the following differences:
% 1) The sample entropy does not calculate self-matching, and the possible log(0) problem is avoided by taking the logarithm in the last step;
% 2) Sample entropy does not depend on the length of the data as does approximate entropy.
% dim: embedding dimension (usually 1 or 2)
% r: Similarity tolerance (usually 0.1*Std(data)~0.25*Std(data))
% data: Time series data. Data must be a 1xN matrix
% tau: Downsampling delay (users can ignore this if the default value is 1)

if nargin < 4, tau = 1;  end
if tau > 1, data = downsample(data, tau);  end

N = length(data);
result = zeros(1,2);

for m = dim:dim+1
Bi = zeros(N-m+1,1);
dataMat = zeros(N-m+1,m);

% Sets the data matrix to be constructed as an m-dimensional vector
for i = 1:N-m+1
dataMat(i,:) = data(1,i:i+m-1);
end

% Calculate the number of similar modes using distance
for j = 1:N-m+1
% Calculates Chebyshev distance, excluding self-matching cases
Dist = Max (abs (dataMat - repmat (dataMat (j:), N - m + 1, 1)), [], 2);
% Statistics the number of dist less than or equal to r
D = (dist <= r);
% does not include self-matching cases
Bi(j,1) = (sum(D)-1)/(N-m);
end

% Find the mean of all Bi's
result(m-dim+1) = sum(Bi)/(N-m+1);

end
% The calculated sample entropy
sampEn = -log(result(2)/result(1));

end