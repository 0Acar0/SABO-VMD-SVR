%Approximate entropy function
function ApEn = kApproximateEntropy(data, dim, r)
% Computational approximate entropy (ApEn)
% data - The data to be analyzed needs to be one-dimensional data
% dim-mode dimension
% r- Threshold size, a boat select r=0.1-0.25,
data = data(:);  % forces data to be converted to column orientation
N = length(data); % Gets the data length
result = zeros(1,2); % Indicates the initialization parameter
for j = 1:2
m = dim+j-1;
C = zeros(1,N-m+1); %C value initialization,
dataMat = zeros (m, N-m+1); % initialization
for i = 1:m
dataMat(i,:) = data(i:N-m+i); %
end
% counting similar patterns using distance calculation
for i = 1:N-m+1
tempMat = abs(dataMat - repmat (dataMat(:,i),1,N-m+1));
boolMat = any((tempMat>r),1);
C(i)=sum(~boolMat)/(N-m+1);
end
phi(j) = sum(log(C))/(N-m+1);
end
ApEn = phi(1)-phi(2);
end
