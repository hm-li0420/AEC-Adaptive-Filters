function [e,y, w] = NLMS(d, x, M, aplha, beta)
    % 输入:
    % d  - 麦克风语音
    % x  - 远端语音
    % a  - 偏置参数
    % M  - 滤波器阶数
    
    % 输出:
    % e - 近端语音估计
    % y - 远端回声估计
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

    xx = zeros(M,1);
    w1 = zeros(M,1);        % 滤波器权重
    y = zeros(d_length,1);   % 近端语音
    e = zeros(d_length,1);  % 误差

    for n = 1:d_length
        %在输入信号x前补上M-1个0，使输出y与输入具有相同长度
        xx = [xx(2:M);x(n)];    % (39,1)+(1,1)=(40,1)
        y(n) = w1' * xx;        % (40,1)'*(40,1)=1
        mu = aplha/(beta + xx'*xx);    % 步长
        e(n) = d(n) - y(n);     % 误差
        
        if e(n) > 2
            e(n)=2;
        elseif e(n)<-2
            e(n)=-2;
        end
        
        w1 = w1 + mu * e(n) * xx;
        w(:,n) = w1;
    end
end
