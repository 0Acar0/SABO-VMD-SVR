%% permutation entropy algorithm
function [pe ,hist] = PermutationEntropy(y,m,t)

%  Calculate the permutation entropy(PE)


%  Input:   y: time series;
% m: order of permuation entropy Embedding dimension
% t: delay time of permuation entropy

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
% πÈ“ªªØ
pe=pe/log(factorial(m));
end
