% 假设你的原始数据存储在一个名为data的向量中，每年12个数据点
% data应该是一个列向量，其中每行对应一个时间步长

% 假设你有一些预测的时间步数，比如外推一年 (12个月)
forecast_steps = 12;

% 搜索最佳差分阶数和季节性阶数
best_aic = Inf;
best_model = [];
for d = 0:2 % 最大差分阶数为2
    for s = 1:12 % 季节性周期为1到12
        model = arima('D', d, 'Seasonality', s);
        try
            model = estimate(model, data);
            if model.ModelCriterion.AIC < best_aic
                best_aic = model.ModelCriterion.AIC;
                best_model = model;
            end
        catch
            continue; % 如果模型拟合失败，继续下一个组合
        end
    end
end

% 使用得到的最佳模型进行预测
forecasted_data = forecast(best_model, forecast_steps, 'Y0', data);

% 显示结果
time_steps = 1:length(data);
future_time_steps = (length(data) + 1):(length(data) + forecast_steps);

figure;
plot(time_steps, data, 'b', 'LineWidth', 1.5);
hold on;
plot(future_time_steps, forecasted_data, 'r--', 'LineWidth', 1.5);
legend('原始数据', '预测数据');
xlabel('时间步长');
ylabel('数据');
title('季节性时间序列预测（自动季节性ARIMA）');
grid on;
