% %
%clc
%clear
fs=12000; % sampling frequency
Ts=1/fs; % sampling period
L=120; % sampling number
t=(0:L-1)*Ts; % time series
STA=1;  % Sampling start position
% -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- import the inner ring fault data -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
%load 105.mat
%X = X105_DE_time(1:L)';  Here you can choose DE(driver end acceleration), FE(fan end acceleration), BA(base acceleration), directly change the variable name, choose one.
% -- -- -- -- -- -- -- -- -- some sample parameters forVMD: for VMD sample parameter Settings -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
alpha = 2500;        % moderate bandwidth constraint: Moderate bandwidth constraint/penalty factor
tau = 0;           % noise-tolerance (no strict fidelity enforcement) : Noise -tolerance (no strict fidelity enforcement)
K = 10;               % modes: The number of decomposed modes, which can be set by yourself. 8 is used as an example.
DC = 0;              % no DC part imposed: No DC part imposed
init = 1;            % initialize omegas uniformly: omegas are uniformly initialized
tol = 1e-7;
% -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- the Run actual VMD code: data VMD decomposition -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
[u, u_hat, omega] = VMD(X, alpha, tau, K, DC, init, tol);  % Where u is the IMF component obtained by decomposition
figure(1);
imfn=u;
n=size(imfn,1);
Subplot (n + 1,1,1);
plot(t,X);  % fault signal
ylabel(' original signal ','fontsize',12,'fontname',' Typeface ');

for n1=1:n
Subplot (n + 1, 1, n1 + 1);
plot(t,u(n1,:)); % Output IMF component, where a(:,n) represents the NTH column element of matrix a, and u(n1,:) represents the n1 row element of matrix u
ylabel(['IMF' int2str(n1)]); %int2str(i) is the rounding of the value i to a character, named on the Y-axis
end
xlabel(' time \itt/s','fontsize',12,'fontname',' Typeface ');
% -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- computing center frequency determine the number of K -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
average=mean(omega);   % averages omega, which is the center frequency.
figure(2)
for i = 1:K
Hy(i,:)= abs(hilbert(u(i,:)));
subplot(K,1,i);
plot(t,u(i,:),'k',t,Hy(i,:),'r');
xlabel(' sample point '); ylabel(' amplitude ')
grid;  legend(' signal ',' envelope ');
end
title(' signal and envelope ');
set(gcf,'color','w');
%% Draw the envelope spectrum
figure('Name',' envelope spectrum ','Color','white');
nfft=fix(L/2);
for i = 1:K
p=abs(fft(Hy(i,:)));  % and fft, get p, is the envelope fft- envelope spectrum
p = p/length(p)*2;
p = p(1: fix(length(p)/2));
subplot(K,1,i);
plot((0: nFFT-1)/nfft*fs/2,p) % Plots the envelope spectrum
%xlim([0 600]) % shows the low band of the envelope spectrum, and this code can choose whether to comment or not according to the situation
if i ==1
title(' envelope spectrum '); xlabel(' Frequency '); ylabel(' amplitude ')
else
xlabel(' Frequency '); ylabel(' amplitude ')
end
end
set(gcf,'color','w');
%% Calculates the kurtosis value
for i=1:K
a(i)=kurtosis(imfn(i,:)); % kurtosis
disp(['IMF',num2str(i),' the kurtosis value is: ',num2str(a(i))])
end
figure
b = bar(a,0.3);
xlabel(' Modal function '); ylabel(' kurtosis value ')
set(gca,'xtick',1:1:K);
set(gca,'xticklabel',{'IMF1','IMF2','IMF3','IMF4','IMF5','IMF6','IMF7','IMF8'});
xtips1 = b.XEndPoints;
ytips1 = b.YEndPoints;  % Gets the XEndPoints and YEndPoints properties of the Bar object
labels1 = string(b.YData);  % Gets the bar end coordinates
text(xtips1, ytips1, labels1, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
% Specifies vertical and horizontal alignment so that the value appears centered at the end of the bar
%% Energy entropy
for i=1:K
Eimf(i) = sum(imfn(i,:).^2,2);
end
disp(['IMF weight of energy '])
disp(Eimf(1:K))
% energy entropy
E = sum(Eimf);
for j = 1:K
p(j) = Eimf(j)/E;
HE(j)=-sum(p(j).*log(p(j)));
end
disp('EMD energy entropy =%.4f');
disp(HE(1:K));



%% Approximate entropy call
for i1=1:K
imf0=imfn(i1,:); %
x=imf0';
ApEnx = kApproximateEntropy(x, 2, 0.15); % approximate entropy calculation
ApEn(1,i1)= ApEnx;
Disp ([' the IMF ', num2str (i1), 'approximate entropy is:', num2str (ApEn (1, i1))))
end


%% envelope entropy
for i = 1:K
xx= abs(hilbert(u(i,:)));  % Minimum envelope entropy calculation formula!
xxx = xx/sum(xx);
ssum=0;
for ii = 1:size(xxx,2)
bb = xxx(1,ii)*log(xxx(1,ii));
ssum=ssum+bb;
end
Enen(i,:) = -ssum; % Envelope entropy for each IMF component
Disp ([' the IMF ', num2str (I), 'the envelope of entropy is:', num2str (Enen (I, :)))
end
ff = min(Enen); % to obtain the local minimum envelope entropy, generally use intelligent optimization algorithm to optimize VMD, the minimum fitness function is used
disp([' local minimum envelope entropy is: ',num2str(ff)])

%% spectral diagram
figure('Name',' spectral diagram ','Color','white');
for i = 1:K
p=abs(fft(u(i,:)));
subplot(K,1,i);
plot((0:L-1)*fs/L,p)
xlim([0 fs/2])
if i ==1
title(' Spectrogram '); xlabel(' Frequency '); ylabel(['IMF' int2str(i)]); %int2str(i) is the rounding of the value i to a character, named on the Y-axis
else
xlabel(' Frequency '); ylabel(['IMF' int2str(i)]); %int2str(i) is the rounding of the value i to a character, named on the Y-axis
end
end
set(gcf,'color','w');


%% fuzzy entropy call
R0 = 0.15; % r is the limit of similarity tolerance
for i1=1:K
imf0=imfn(i1,:); %
x=imf0;
r=r0*std(x);
FuEnx = FuzzyEntropy(x, 6, r,2,1); % fuzzy entropy
FuEn(1,i1)= FuEnx;
Disp ([' the IMF ', num2str (i1), 'fuzzy entropy is:', num2str (FuEn (1, i1))))
end

%% permutation entropy call
M = 3;   % Embedding dimension
T = 1;   % delay time
for i1=1:K
imf0=imfn(i1,:); %
x=imf0;
PerEnx = PermutationEntropy(x,M,T);     % permutation entropy calculation
PerEn(1,i1)= PerEnx;
Disp ([' the IMF ', num2str (i1), 'permutation entropy is:', num2str (PerEn (1, i1))))
end


%% Multi-scale entropy call
m = 6;           %m: order of permuation entropy;
t = 1;           %t: delay time of permuation entropy;
scale = 30;      %scale: the scale factor;
for i1=1:K
imf0=imfn(i1,:); %
x=imf0';
MPE = MultiscalePermutationEntropy(x,m,t,scale);     % Multi-scale entropy calculation
B(i1,:)= MPE;
end
disp(' Multi-scale arrangement entropy is: ')
disp(num2str(B))

%% Sample entropy call
dim = 2;    % dim: embedding dimension (usually 1 or 2)
tau = 1;    % Downsampling delay time (Users can ignore this if the default value is 1)
for i1=1:K
imf0=imfn(i1,:); %
x=imf0;
r = 0.2*std(x); % r: Similarity tolerance (usually 0.1*Std(data)~0.25*Std(data))
SaEnx = SampleEntropy(dim,r, x, tau);     % sample entropy calculation
SaEn(1,i1)= SaEnx;
Disp ([' the IMF ', num2str (i1), 'the sample entropy is:', num2str (SaEn (1, i1))))
end