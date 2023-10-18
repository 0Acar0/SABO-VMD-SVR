clear
clc
close all
xz = 1;  %xz=1 or 2, Select 1, taking the minimum envelope entropy as the fitness function; 
              select 2, taking the minimum sample entropy as the fitness function; 
              select 3, taking the minimum information entropy as the fitness function; 
              select 4, taking the minimum permutation entropy as the fitness function
if xz == 1  
    fobj=@EnvelopeEntropyCost;          %Minimum envelope entropy
elseif xz == 2
    fobj=@SampleEntropyCost;            %Minimum sample entropy
elseif xz == 3
    fobj=@infoEntropyCost;              %Minimum information entropy
elseif xz == 4
    fobj=@PermutationEntropyCost;       %Minimum permutation entropy
end

%% read data
%load 105.mat
%da = X105_DE_time(6001:7000); %The data format is n rows x 1 column. The number of columns must be 1.

%% Set parameters
lb = [100 3];    %Lower bound for the penalty factor and K
ub = [2500 10];  %The penalty factor and the upper limit of K
dim = 2;            % optimize  Number of  variables
Max_iter=20;       % Maximum number of iterations
SearchAgents_no=30;       %Population size

%% Call the SABO function
[fMin , bestX, Convergence_curve ] = SABO(SearchAgents_no,Max_iter,lb,ub,dim,fobj,da);

%% Draw a fitness function graph and output the best parameters
figure
plot(Convergence_curve,'Color',[0.9 0.5 0.1],'Marker','>','LineStyle','--','linewidth',1);
title('Objective space')
xlabel('Iteration');
ylabel('Best score obtained so far');
legend('SABO”≈ªØVMD')
display(['The best solution obtained by SABO is : ', num2str(round(bestX))]);  %Output optimum position
display(['The best optimal value of the objective funciton found by SABO is : ', num2str(fMin)]);  %Output the optimum fitness value



