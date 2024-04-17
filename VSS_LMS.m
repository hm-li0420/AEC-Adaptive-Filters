% 学习率的上下界
% 以上信息可从论文中获取，这些值是在下述论文中定义的
% R. H. Kwong and E. W. Johnston, "A variable step size LMS algorithm," in IEEE Transactions on Signal Processing, vol. 40, no. 7, pp. 1633-1642, July 1992.
function [e, y, w] = VSS_LMS(d, x, M)
% Inputs:
% d  - 麦克风语音
% x  - 远端语音
% mu - 步长，0.05
% M  - 滤波器阶数，也称为抽头数
%
% Outputs:
% e - 输出误差
% y  - 输出系数
% w - 滤波器参数

    d_length = length(d);
    if (d_length <= M)  
        print('error: 信号长度小于滤波器阶数！');
        return; 
    end
    if (d_length ~= length(x))  
        print('error: 输入信号和参考信号长度不同！');
        return; 
    end
    
    % 定义VSS-LMS算法的初始参数
    
    input_var = var(x);     % 输入的方差
    % 如果mu_max和mu_min之间的差异不够大，LMS和VSS LMS算法的误差曲线都是相同的
    mu_max = 1/(input_var*M); % 上界=1/(filter_length * input variance)
    mu_min = 0.0004;        % 下界=LMS算法的学习率
    mu_VSS(1)=1;    % VSS算法的mu初始值
    alpha  = 0.97;
    gamma = 4.8e-4;
    
    xx = zeros(M,1);
    w1 = zeros(M,1);    % 滤波器权重
    y = zeros(d_length,1);  % 近端语音
    e = zeros(d_length,1);  % 误差

    for n = 1:d_length
        xx = [xx(2:M);x(n)];    % 纵向拼接 或者[x(n); xx(1:end-1)]
        y(n) = w1' * xx;        % (40,1)'*(40,1)=1; (73113,1)
        e(n) = d(n) - y(n);
        w1 = w1 + mu_VSS(n) * e(n) * xx;   % 更新权重系数 (40,1)
        w(:,n) = w1;        % (40, 73113)
        
        mu_VSS(n+1) = alpha * mu_VSS(n) + gamma * e(n) * e(n) ;% 使用VSS算法更新mu值
       % 检查论文中给出的mu的约束条件
        if (mu_VSS(n+1)>mu_max)
            mu_VSS(n+1)=mu_max; % max
        elseif(mu_VSS(n+1)<mu_min)
            mu_VSS(n+1)= mu_min;
        else
            mu_VSS(n+1) = mu_VSS(n+1) ;
        end
    end
    % 和上面类似
%     for n = M:d_length
%         xx = x(n:-1:n-M+1);    % 纵向拼接  (40~1)-->(41~2)-->(42~3)....
%         y(n) = w1' * xx;        % (40,1)'*(40,1)=1; (73113,1)
%         e(n) = d(n) - y(n);
%         w1 = w1 + mu * e(n) * xx;   % (40,1)
%         w(:,n) = w1;        % (40, 73113)
%     end
end