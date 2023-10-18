% �������ԭʼ���ݴ洢��һ����Ϊdata�������У�ÿ��12�����ݵ�
% dataӦ����һ��������������ÿ�ж�Ӧһ��ʱ�䲽��

% ��������һЩԤ���ʱ�䲽������������һ�� (12����)
forecast_steps = 12;

% ������Ѳ�ֽ����ͼ����Խ���
best_aic = Inf;
best_model = [];
for d = 0:2 % ����ֽ���Ϊ2
    for s = 1:12 % ����������Ϊ1��12
        model = arima('D', d, 'Seasonality', s);
        try
            model = estimate(model, data);
            if model.ModelCriterion.AIC < best_aic
                best_aic = model.ModelCriterion.AIC;
                best_model = model;
            end
        catch
            continue; % ���ģ�����ʧ�ܣ�������һ�����
        end
    end
end

% ʹ�õõ������ģ�ͽ���Ԥ��
forecasted_data = forecast(best_model, forecast_steps, 'Y0', data);

% ��ʾ���
time_steps = 1:length(data);
future_time_steps = (length(data) + 1):(length(data) + forecast_steps);

figure;
plot(time_steps, data, 'b', 'LineWidth', 1.5);
hold on;
plot(future_time_steps, forecasted_data, 'r--', 'LineWidth', 1.5);
legend('ԭʼ����', 'Ԥ������');
xlabel('ʱ�䲽��');
ylabel('����');
title('������ʱ������Ԥ�⣨�Զ�������ARIMA��');
grid on;
