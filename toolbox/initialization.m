function Positions = initialization(SearchAgents_no, dim, ub, lb)

%% 优化参数个数
Boundary_no= size(ub, 2);

%% 只有一个优化参数
if Boundary_no == 1
    Positions = rand(SearchAgents_no, dim) .* (ub - lb) + lb;
end

%% 多个优化参数，且上下界不同
if Boundary_no > 1
    for i = 1:dim
        ub_i = ub(i);
        lb_i = lb(i);
        Positions(:, i) = rand(SearchAgents_no, 1) .* (ub_i - lb_i) + lb_i;
    end
end