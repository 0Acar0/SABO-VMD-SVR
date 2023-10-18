
function[Best_score,Best_pos,SABO_curve]=SABO(N,T,lo,hi,m,fitness,da)

lo=ones(1,m).*(lo);                              % Lower limit for variables
hi=ones(1,m).*(hi);                              % Upper limit for variables

%% INITIALIZATION
for i=1:m
    X(:,i) = lo(i)+rand(N,1).*(hi(i) - lo(i));                          % Initial population
end

for i =1:N
    L=X(i,:);
    fit(i)=fitness(L,da);
end
%%

for t=1:T  % algorithm iteration
    
    %%  update: BEST proposed solution
    [Fbest , blocation]=min(fit);
    if t==1
        xbest=X(blocation,:);                                           % Optimal location
        fbest=Fbest;                                           % The optimization objective function
    elseif Fbest<fbest
        fbest=Fbest;
        xbest=X(blocation,:);
    end
    %%
    DX=zeros(N,m);  
    for i=1:N
        %% based om Eq(4)
        for j=1:N
            I=round(1+rand+rand);
            for d=1:m
                DX(i,d)=DX(i,d)+(X(j,d)-I.*X(i,d)).*sign(fit(i)-fit(j));
            end
        end
        X_new_P1= X(i,:)+((rand(1,m).*DX(i,:))./(N));
        X_new_P1 = max(X_new_P1,lo);
        X_new_P1 = min(X_new_P1,hi);
        
        %% update position based on Eq (5)
        L=X_new_P1;
        fit_new_P1=fitness(L,da);
        if fit_new_P1<fit(i)
            X(i,:) = X_new_P1;
            fit(i) = fit_new_P1;
        end
        %%            
    end% end for i=1:N
    %%
    
    best_so_far(t)=fbest;
    average(t) = mean (fit);
    
    disp(['SABO: At iteration ', num2str(t), ' ,the best fitness is ', num2str(fbest)])
    disp(['SABO: At iteration ',num2str(t),'the best position is£º[',num2str(fix(xbest)),']'])
end
Best_score=fbest;
Best_pos=xbest;
SABO_curve=best_so_far;
end

