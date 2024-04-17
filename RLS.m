function [e, y, w] = RLS(d, x, M, lamda)
% Inputs:
% d  - 麦克风语音
% x  - 远端语音
% lamda - the weight parameter, 权重
% M  - the number of taps. 滤波器阶数
%
% Outputs:
% e - 大小为Ns的输出误差向量
% y - 输出的近端语音
% w - 滤波器参数

    Ns = length(d);
    if (Ns <= M)  
        print('error: 信号长度小于滤波器阶数！');
        return; 
    end
    if (Ns ~= length(x))
        print('error: 输入信号和参考信号长度不同！');
        return;
    end

    I = eye(M);
    a = 0.01;
    p = a * I;
    
    xx = a.*zeros(M,1);
    w1 = zeros(M,1);        % 滤波器权重
    y = zeros(Ns, 1);        % 近端语音
    e = zeros(Ns, 1);       % 误差
    
    for n = 1:Ns
        %在输入信号x后补上M-1个0，使输出y与输入具有相同长度
        xx = [x(n); xx(1:M-1)];
        k = (p * xx) ./ (lamda + xx' * p * xx + 1e-9);
        
        if k > 2000
            k=2000;
        elseif k<-2000
            k=-2000;
        end
        
        y(n) = xx'*w1;
        e(n) = d(n) - y(n);
        
        if e(n) > 2
            e(n)=2;
        elseif e(n)<-2
            e(n)=-2;
        end
        
        w1 = w1 + k * e(n);
        p = (p - k * xx' * p) ./ (lamda+ 1e-9);
        
        if p > 10
            p=10;
        elseif p<-20
            p=-20;
        end
        
        w(:,n) = w1;
        
        if w1 > 2
            w1=2;
        elseif w1<-2
            w1=-2;
        end
        
    end
end
