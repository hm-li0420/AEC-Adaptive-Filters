clear; clc;

%% Get far-end and near-end speech files
% near-end
[filename1, pathname1] = uigetfile('*.wav', 'Select near-end speech file');
Filename1 = strcat(pathname1, filename1);
[v, fs1] = audioread(Filename1);
% far-end
[filename2, pathname2] = uigetfile('*.wav', 'Select far-end speech file');
Filename2 = strcat(pathname2, filename2);
[x, fs2] = audioread(Filename2);

%% Create echo in a room
M = fs1/2 + 1;
% creating impuls response
[B,A] = cheby2(4,20,[0.1 0.7]);
impulseResponseGenerator = dsp.IIRFilter('Numerator', [zeros(1,6) B], 'Denominator', A);
% creating a room
roomImpulseResponse = impulseResponseGenerator((log(0.99*rand(1,M)+0.01).*sign(randn(1,M)).*exp(-0.002*(1:M)))');
roomImpulseResponse = roomImpulseResponse/norm(roomImpulseResponse)*4;
room = dsp.FIRFilter('Numerator', roomImpulseResponse');
% far-end echo speech
x_echo = room(x);
r = x + x_echo*0.3;
% % plot
% figure;
% subplot(3,1,1); plot(x); title('far-end speech');
% subplot(3,1,2); plot(x_echo); title('echo');
% subplot(3,1,3); plot(r); title('far-end echo speech');

%% Combine far-end echo speech and near-end speech
d = r + v;

%% Double-talk detector -- Geigel Algorithm
% look up time, N in Equation (3.23)
N = fs1/50;
% threshold, T in Equation (3.23)
% need to interactively change according to different speech 
T = 0.5;
dtd = zeros(1,length(d));
% Geigel Algorithm over the whole signal
for index = 1:length(d)
    if index-N+1 >= 1
        sub = abs(x((index-N+1):index));
    else
        sub = abs(x(1):x(index));
    end
    Geigel = max(sub)/abs(d(index));
    if Geigel < T
        dtd(index) = 1;    % presence of near-end speech
    else
        dtd(index) = 0;    % absent of near-end speech
    end
end

%% NLMS adaptive filter
% initial condition
mu = 0.00001;
gamma = 10;
w = zeros(1,length(d));
y = zeros(1,length(d));
e = zeros(1,length(d));
xT = transpose(x);
% loop
for index = 2:length(d)
    if dtd(index) == 0
        y(index) = xT(index) .* w(index);
        e(index) = d(index) - y(index);
        w(index+1) = w(index) + e(index).*x(index).*2.*mu./(gamma + x(index).*xT(index));
    end
end
soundsc(y,fs1);

% %% Non-linear processor
% threshold = 0.00001;   % need to change interactively according to different speech 
% result = zeros(1,length(y));
% for index = 1:length(y)
%     if abs(y(index)) < threshold
%         result(index) = y(index)./10;
%     else
%         result(index) = y(index);
%     end
% end
% figure;
% subplot(2,1,1); plot(x); title('original far-end speech');
% subplot(2,1,2); plot(result); title('echo cancellation result');
% 
% %% Echo Return Loss Enhancement (ERLE)
% error_signal = abs(result-x');
% power_d = sum(d.^2);
% power_error_signal = sum(error_signal.^2);
% ERLE = 10*log10(power_d./power_error_signal);