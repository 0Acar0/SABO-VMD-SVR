function MPE = MultiscalePermutationEntropy(data,m,t,scale)
%  Calculate the Multiscale Permutation Entropy (MPE)
%  Input:   data: time series;
%           m: order of permuation entropy;
%           t: delay time of permuation entropy; 
%           scale: the scale factor;
%  Output: 
%           MPE: multiscale permuation entropy.

MPE=[];

for j=1:scale
    Xs = Multi(data,j);
    PE = PermutationEntropy(Xs,m,t);
    MPE=[MPE,PE];
end

end

function M_Data = Multi(data,scale)
%  generate the consecutive coarse-grained time series
%  Input:   data: time series;
%           scale: the scale factor
%  Output: 
%           M_Data: the coarse-grained time series at the scale factor S
L = length(data);
J = fix(L/scale);

for i=1:J
    M_Data(i) = mean(data((i-1)*scale+1:i*scale));
end
end

function [apen] = PermutationEntropy(data,m,t)
%  Calculate the permutation entropy
%  Input:   
%           data: time series;
%           m: order of permuation entropy
%           t: delay time of permuation entropy, 
%  Output:  
%           apen: permuation entropy
%Ref: G Ouyang, J Li, X Liu, X Li, Dynamic Characteristics of Absence EEG Recordings with Multiscale Permutation %  
%                             Entropy Analysis, Epilepsy Research, doi: 10.1016/j.eplepsyres.2012.11.003
%     X Li, G Ouyang, D Richards, Predictability analysis of absence seizures with permutation entropy, Epilepsy %  
%                            Research,  Vol. 77pp. 70-74, 2007
N = length(data);
permlist = perms(1:m);
c(1:length(permlist))=0;
    
 for i=1:N-t*(m-1)
     [~,iv]=sort(data(i:t:i+t*(m-1)));
     for j=1:length(permlist)
         if (permlist(j,:)-iv)==0 %如果if后的判定表达式为(逻辑)向量，则判断其中所有元素是否全非0
             c(j) = c(j) + 1 ;
         end
     end
 end

hist = c; 
c = hist(hist~=0);
p = c/sum(c);
pe = -sum(p .* log(p));
% normalized
apen=pe/log(factorial(m));
end
