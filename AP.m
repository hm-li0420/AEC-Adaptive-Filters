% 参考：
% https://github.com/rohillarohan/System-Identification-for-Echo-Cancellation
% https://github.com/3arbouch/ActiveNoiseCancelling/blob/c29ccfcd869657a5f58e1ce7755fe08c3a195bd9/ANC/BookExamples/MATLAB/APLMS_AEC_mono.m
function [e,y,w] = AP(d,x,mu,M,L,psi)
    % input: 
    % d   -- 麦克风信号
    % x   -- 远端语音
    % mu  -- 步长
    % M   -- 滤波器阶数 40
    % L   -- 列数 2
    % psi -- 1e-4
    %
    % Outputs:
    % e - 输出误差 the output error vector of size Ns
    % y  - 输出系数 output coefficients
    % w - 滤波器参数 filter parameters
    
    d_length = length(d);
    if (d_length <= M)  
        print('error: 信号长度小于滤波器阶数！');
        return; 
    end
    if (d_length ~= length(x))  
        print('error: 输入信号和参考信号长度不同！');
        return; 
    end
    
    XAF = zeros(M,L);
    w1 = zeros(M,1);  % 滤波器权重 (40, 1)
    y = zeros(d_length,1);  % 估计的近端语音
    e = zeros(d_length,1);  % 误差
    
    for m = M+L:d_length    % 采样点数
        for k = 1:L % 列数
            XAF(:,k) = x(m-k+1:-1:m-k+1-M+1);   % (40,2)
        end
%         y(m) = (XAF'*w)'*(XAF'*w);        % 不太确定是不是这样(40,2)'*(40,1)=(2,1)
        E = d(m:-1:m-L+1) -XAF'*w1;    % (2,1)-(2,1)
        w1 = w1 + mu*XAF*inv((XAF'*XAF + psi*eye(L)))*E;
        w(:,m) = w1;        % (40, 73113)
        e(m) = E(1)'*E(1);
        
        
        if e(m) > 2
            e(m)=2;
        elseif e(m)<-2
            e(m)=-2;
        end
    end
end