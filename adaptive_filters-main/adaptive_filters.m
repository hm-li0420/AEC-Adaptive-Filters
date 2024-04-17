%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameter Initialization
nMCT = 3;                         %Number of Monte-Carlo Trials
M = 2;                            %adaptive filter length
N = 5000;                         %number of samples
sigma_V_s = 1;                    %standard deviation of V_s[n]
sigma_V_i = 0.001;                  %standard deviation of V_i[n]
sigma_V_e = 0;                  %standard deviation of V_e[n]
H_d = [500 30];                      % H_d
a_s = 0.6;                          %filter coefficient of H_s
 
%% Step Sizes
mu_NLMS = 0.01;                     %step size of NLMS
mu_PNLMS = 0.01;                     %step size of PNLMS
mu_IPNLMS = 0.01;                     %step size of IPNLMS
mu_GCIPNLMS = 0.01;                     %step size of GCIPNLMS
 
%% NLMS Parameters
Epsilon_NLMS = 0.01;        %Epsilon parameter for NLMS algorithm
 
%% PNLMS Parameters
Epsilon_PNLMS = 0.01;
g_PNLMS = zeros(1,M);
g_bar_PNLMS = 1; 
sigma_PNLMS = 0.01;
ro_PNLMS = 0.01;
 
%% IPNLMS Parameters
K_IPNLMS = eye(M);
g_IPNLMS = zeros(1,M);
g_bar_IPNLMS = 1; 
sigma_IPNLMS = 0.01;
ro_IPNLMS = 0.01;
alpha_IPNLMS = -0.5; %between -1 and 1
% Epsilon_IPNLMS = 0.01;
Epsilon_IPNLMS = ((1-alpha_IPNLMS)/(2*M))*Epsilon_NLMS;
eps_IPNLMS = 0.0001;
 
%% GCIPNLMS Parameters
Q_GCIPNLMS = eye(M);
g_GCIPNLMS = zeros(1,M);
sigma_GCIPNLMS = 0.01;
alpha_GCIPNLMS = -0.5; %between -1 and 1
beta_GCIPNLMS = 0.9999; %taken from paper
Epsilon_GCIPNLMS = 0.01;
 
%% All filter coefficients for plotting
 
AllFilterCoefficients_NLMS = zeros(nMCT, N-M, M);      %Holds all filter coefficients of NLMS for plotting
 
AllFilterCoefficients_PNLMS = zeros(nMCT, N-M, M);     %Holds all filter coefficients of PNLMS for plotting
AllFilterCoefficients_IPNLMS = zeros(nMCT, N-M, M);     %Holds all filter coefficients of IPNLMS for plotting
AllFilterCoefficients_GCIPNLMS = zeros(nMCT, N-M, M);     %Holds all filter coefficients of GCIPNLMS for plotting
%%    Monte Carlo Trials
 
for i = 1: nMCT 
    h_initial = i*ones(1,M);              %start with different initial conditions
    
    %% Initialize filter coefficients
    h_adapt_NLMS = zeros(N-M,M);     % adaptive filter coefficients for NLMS
    h_adapt_PNLMS = zeros(N-M,M);    % adaptive filter coefficients for PNLMS
    h_adapt_IPNLMS = zeros(N-M,M);    % adaptive filter coefficients for IPNLMS
    h_adapt_GCIPNLMS = zeros(N-M,M);    % adaptive filter coefficients for GCIPNLMS
     
    h_adapt_NLMS(1,:) = h_initial;
 
    %new algorithms for final
    h_adapt_PNLMS = h_initial;
    h_adapt_IPNLMS = h_initial;    
    h_adapt_GCIPNLMS = h_initial;
    %%
    V_s = sigma_V_s*randn(N,1);     %V_s[n]    
    s = filter(1,[1 -a_s], V_s);    %obtain s[n]
    x = s + sigma_V_i*randn(N,1);   %x[n]
    d = filter(H_d,1,s);            %d[n]
    
    %% Initialize Errors
    e_NLMS = zeros(1,N);            %e[n] of NLMS
    e_PNLMS = zeros(1,N);           %e[n] of PNLMS
    e_IPNLMS = zeros(1,N);          %e[n] of IPNLMS
    e_GCIPNLMS = zeros(1,N);        %e[n] of GCIPNLMS
    %% 
    AlgorithmInput = zeros(1,N);    %e[n] + V_e[n]
    V_e = sigma_V_e*randn(N,1);     %V_e[n]
 
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   NLMS Algorithm
    for n = 1:N-M
        y = h_adapt_NLMS(n,:)*x(n:n+M-1); % filter coefficients are reversed
        e_NLMS(n) = d(n+M-1) - y;
        AllFilterCoefficients_NLMS(i, n, :) = h_adapt_NLMS(n,:);
        AlgorithmInput(n) = V_e(n) + e_NLMS(n);
        h_adapt_NLMS(n+1,:) = h_adapt_NLMS(n,:) + (mu_NLMS /(Epsilon_NLMS + norm(x(n:n+M-1))^2))*AlgorithmInput(n)*x(n:n+M-1)';
    end
    h_adapt_final_NLMS = h_adapt_NLMS(end,:);
    h_adapt_final_NLMS = fliplr(h_adapt_final_NLMS);
 
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   PNLMS Algorithm
    for n = 1:N-M
        y = h_adapt_PNLMS(n,:)*x(n:n+M-1); % filter coefficients are reversed
        e_PNLMS(n) = d(n+M-1) - y;
        AllFilterCoefficients_PNLMS(i, n, :) = h_adapt_PNLMS(n,:);
        %Calculation of g and g_bar
        l_inf = max(abs(h_adapt_PNLMS(n,:)));
        l_prime_inf = max(sigma_PNLMS, l_inf);
        
        for gg = 1:M
            g_PNLMS(gg) = max(ro_PNLMS*l_prime_inf, abs(h_adapt_PNLMS(n,gg)));
        end
        
        g_bar_PNLMS = (1/M)*sum(g_PNLMS);
        AlgorithmInput(n) = V_e(n) + e_PNLMS(n);
        
        h_adapt_PNLMS(n+1,:) = h_adapt_PNLMS(n,:) + (flip(g_PNLMS)/g_bar_PNLMS).*((mu_PNLMS /(Epsilon_PNLMS + norm(x(n:n+M-1))^2))*AlgorithmInput(n)*x(n:n+M-1)');
    end
    h_adapt_final_PNLMS = h_adapt_PNLMS(end,:);
    h_adapt_final_PNLMS = fliplr(h_adapt_final_PNLMS);
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   IPNLMS Algorithm
    for n = 1:N-M
        y = h_adapt_IPNLMS(n,:)*x(n:n+M-1); % filter coefficients are reversed
        e_IPNLMS(n) = d(n+M-1) - y;
        AllFilterCoefficients_IPNLMS(i, n, :) = h_adapt_IPNLMS(n,:);
        
        for gg = 1:M
            K_IPNLMS(gg,gg) = (1 - alpha_IPNLMS)/(2*M) + (1 + alpha_IPNLMS)*(abs(h_adapt_IPNLMS(n,gg))/(2 * sum(h_adapt_IPNLMS(n,:)) + eps_IPNLMS));
        end
        
        AlgorithmInput(n) = V_e(n) + e_IPNLMS(n);
        h_adapt_IPNLMS(n+1,:) = h_adapt_IPNLMS(n,:) + ((mu_IPNLMS /(Epsilon_IPNLMS + x(n:n+M-1)'*K_IPNLMS*x(n:n+M-1)))*AlgorithmInput(n)*x(n:n+M-1)'*K_IPNLMS);
    end
    h_adapt_final_IPNLMS = h_adapt_IPNLMS(end,:);
    h_adapt_final_IPNLMS = fliplr(h_adapt_final_IPNLMS);
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   GCIPNLMS Algorithm
    for n = 1:N-M
        y = h_adapt_GCIPNLMS(n,:)*x(n:n+M-1); % filter coefficients are reversed
        e_GCIPNLMS(n) = d(n+M-1) - y;
        AllFilterCoefficients_GCIPNLMS(i, n, :) = h_adapt_GCIPNLMS(n,:);
        
        g_GCIPNLMS = beta_GCIPNLMS*g_GCIPNLMS - (1 - beta_GCIPNLMS)*(x(n:n+M-1)')*inv(x(n:n+M-1)'*x(n:n+M-1))*e_GCIPNLMS(n);
            
        for gg = 1:M
            Q_GCIPNLMS(gg,gg) = (1 - alpha_GCIPNLMS)/(2*M) + (1 + alpha_GCIPNLMS)*(abs(g_GCIPNLMS(gg))/(2 * sum(g_GCIPNLMS) + Epsilon_GCIPNLMS));
        end
        
        AlgorithmInput(n) = V_e(n) + e_GCIPNLMS(n);
        h_adapt_GCIPNLMS(n+1,:) = h_adapt_GCIPNLMS(n,:) + ((mu_GCIPNLMS /(Epsilon_GCIPNLMS + x(n:n+M-1)'*Q_GCIPNLMS*x(n:n+M-1)))*AlgorithmInput(n)*x(n:n+M-1)'*Q_GCIPNLMS);
    end
    h_adapt_final_GCIPNLMS = h_adapt_GCIPNLMS(end,:);
    h_adapt_final_GCIPNLMS = fliplr(h_adapt_final_GCIPNLMS);
  
end
 
%% NLMS ERROR CALCULATION
figure
plot(1:N,abs(e_NLMS).^2,'LineWidth',2)
title('Error Power of NLMS')
xlabel('n')
ylabel('e[n]^2')
grid on
hold on
 
%%Convergence Time Calculation
ConvergenceIndex_NLMS = 0;
ConvergenceBound_NLMS = 10;
for i = numel(e_NLMS):-1:2
    if(abs(e_NLMS(i))^2 > ConvergenceBound_NLMS)
         plot(i+1,abs(e_NLMS(i+1)).^2,'o','LineWidth',4);
         ConvergenceIndex_NLMS = i+1;
         break  
    end
end
%Print convergence index
hold off;
disp(['ConvergenceIndex of NLMS = ',num2str(ConvergenceIndex_NLMS)]);
 
%Convergence Time Indication
figure
for i = 1:nMCT
    plot(AllFilterCoefficients_NLMS(i,:,2),AllFilterCoefficients_NLMS(i,:,1),'-o')
    hold on   
end
title('Convergence of Filter Coefficients for NLMS')
grid on
hold off;
 
%Print other results
% fprintf('Optimal Coefficients of NLMS = %3s\n',h_adapt_final_NLMS);
SteadyStateError_NLMS = mean(e_NLMS(numel(e_NLMS)-20:end));
disp(['Steady State Error of NLMS = ',num2str(SteadyStateError_NLMS)]);
MinimumError_NLMS = min(e_NLMS);
disp(['Min Error of NLMS = ',num2str(MinimumError_NLMS)]);
Error_xs_NLMS = (SteadyStateError_NLMS - MinimumError_NLMS);
disp(['Error_xs = ',num2str(Error_xs_NLMS)]);
Misadjustment_NLMS = Error_xs_NLMS/MinimumError_NLMS;
disp(['Misadjustment of NLMS = ',num2str(Misadjustment_NLMS)]);
 
 
%% PNLMS ERROR CALCULATION
figure
plot(1:N,abs(e_PNLMS).^2,'LineWidth',2)
title('Error Power of PNLMS')
xlabel('n')
ylabel('e[n]^2')
grid on
hold on
 
%%Convergence Time Calculation
ConvergenceIndex_PNLMS = 0;
ConvergenceBound_PNLMS = 10;
for i = numel(e_PNLMS):-1:2
    if(abs(e_PNLMS(i))^2 > ConvergenceBound_PNLMS)
         plot(i+1,abs(e_PNLMS(i+1)).^2,'o','LineWidth',4);
         ConvergenceIndex_PNLMS = i+1;
         break  
    end
end
%Print convergence index
hold off;
disp(['ConvergenceIndex of PNLMS = ',num2str(ConvergenceIndex_PNLMS)]);
 
%Convergence Time Indication
figure
for i = 1:nMCT
    plot(AllFilterCoefficients_PNLMS(i,:,2),AllFilterCoefficients_PNLMS(i,:,1),'-o')
    hold on    
end
title('Convergence of Filter Coefficients of PNLMS')
grid on
hold off;
 
%Print other results
% fprintf('Optimal Coefficients of SignedLMS = %3s\n',h_adapt_final_SignedLMS);
SteadyStateError_PNLMS = mean(e_PNLMS(numel(e_PNLMS)-20:end));
disp(['Steady State Error of PNLMS = ',num2str(SteadyStateError_PNLMS)]);
MinimumError_PNLMS = min(e_PNLMS);
disp(['Min Error of PNLMS = ',num2str(MinimumError_PNLMS)]);
Error_xs_PNLMS = (SteadyStateError_PNLMS - MinimumError_PNLMS);
disp(['Error_xs = ',num2str(Error_xs_PNLMS)]);
Misadjustment_PNLMS = Error_xs_PNLMS/MinimumError_PNLMS;
disp(['Misadjustment of PNLMS = ',num2str(Misadjustment_PNLMS)]);
 
%% IPNLMS ERROR CALCULATION
figure
plot(1:N,abs(e_IPNLMS).^2,'LineWidth',2)
title('Error Power of IPNLMS')
xlabel('n')
ylabel('e[n]^2')
grid on
hold on
 
%%Convergence Time Calculation
ConvergenceIndex_IPNLMS = 0;
ConvergenceBound_IPNLMS = 10;
for i = numel(e_IPNLMS):-1:2
    if(abs(e_IPNLMS(i))^2 > ConvergenceBound_IPNLMS)
         plot(i+1,abs(e_IPNLMS(i+1)).^2,'o','LineWidth',4);
         ConvergenceIndex_IPNLMS = i+1;
         break  
    end
end
%Print convergence index
hold off;
disp(['ConvergenceIndex of IPNLMS = ',num2str(ConvergenceIndex_IPNLMS)]);
 
%Convergence Time Indication
figure
for i = 1:nMCT
    plot(AllFilterCoefficients_IPNLMS(i,:,2),AllFilterCoefficients_IPNLMS(i,:,1),'-o')
    hold on    
end
title('Convergence of Filter Coefficients of IPNLMS')
grid on
hold off;
 
%Print other results
% fprintf('Optimal Coefficients of SignedLMS = %3s\n',h_adapt_final_SignedLMS);
SteadyStateError_IPNLMS = mean(e_IPNLMS(numel(e_IPNLMS)-20:end));
disp(['Steady State Error of IPNLMS = ',num2str(SteadyStateError_IPNLMS)]);
MinimumError_IPNLMS = min(e_IPNLMS);
disp(['Min Error of IPNLMS = ',num2str(MinimumError_IPNLMS)]);
Error_xs_IPNLMS = (SteadyStateError_IPNLMS - MinimumError_IPNLMS);
disp(['Error_xs = ',num2str(Error_xs_IPNLMS)]);
Misadjustment_IPNLMS = Error_xs_IPNLMS/MinimumError_IPNLMS;
disp(['Misadjustment of IPNLMS = ',num2str(Misadjustment_IPNLMS)]);
 
%% GCIPNLMS ERROR CALCULATION
figure
plot(1:N,abs(e_GCIPNLMS).^2,'LineWidth',2)
title('Error Power of GCIPNLMS')
xlabel('n')
ylabel('e[n]^2')
grid on
hold on
 
%%Convergence Time Calculation
ConvergenceIndex_GCIPNLMS = 0;
ConvergenceBound_GCIPNLMS = 10;
for i = numel(e_GCIPNLMS):-1:2
    if(abs(e_GCIPNLMS(i))^2 > ConvergenceBound_GCIPNLMS)
         plot(i+1,abs(e_GCIPNLMS(i+1)).^2,'o','LineWidth',4);
         ConvergenceIndex_GCIPNLMS = i+1;
         break  
    end
end
%Print convergence index
hold off;
disp(['ConvergenceIndex of GCIPNLMS = ',num2str(ConvergenceIndex_GCIPNLMS)]);
 
%Convergence Time Indication
figure
for i = 1:nMCT
    plot(AllFilterCoefficients_GCIPNLMS(i,:,2),AllFilterCoefficients_GCIPNLMS(i,:,1),'-o')
    hold on    
end
title('Convergence of Filter Coefficients of GCIPNLMS')
grid on
hold off;
 
%Print other results
% fprintf('Optimal Coefficients of SignedLMS = %3s\n',h_adapt_final_SignedLMS);
SteadyStateError_GCIPNLMS = mean(e_GCIPNLMS(numel(e_GCIPNLMS)-20:end));
disp(['Steady State Error of GCIPNLMS = ',num2str(SteadyStateError_GCIPNLMS)]);
MinimumError_GCIPNLMS = min(e_GCIPNLMS);
disp(['Min Error of GCIPNLMS = ',num2str(MinimumError_GCIPNLMS)]);
Error_xs_GCIPNLMS = (SteadyStateError_GCIPNLMS - MinimumError_GCIPNLMS);
disp(['Error_xs = ',num2str(Error_xs_GCIPNLMS)]);
Misadjustment_GCIPNLMS = Error_xs_GCIPNLMS/MinimumError_GCIPNLMS;
disp(['Misadjustment of GCIPNLMS = ',num2str(Misadjustment_GCIPNLMS)]);
 
%% Final Plot to see convergence of all
 
nMCTforPlot = 1;
figure
plot(AllFilterCoefficients_NLMS(nMCTforPlot,:,2),AllFilterCoefficients_NLMS(nMCTforPlot,:,1),'-o')
hold on
plot(AllFilterCoefficients_PNLMS(nMCTforPlot,:,2),AllFilterCoefficients_PNLMS(nMCTforPlot,:,1),'-o')
plot(AllFilterCoefficients_IPNLMS(nMCTforPlot,:,2),AllFilterCoefficients_IPNLMS(nMCTforPlot,:,1),'-o')
plot(AllFilterCoefficients_GCIPNLMS(nMCTforPlot,:,2),AllFilterCoefficients_GCIPNLMS(nMCTforPlot,:,1),'-o')
title('Convergence of All Algorithms')
grid on
hold off;
legend('NLMS','PNLMS','IPNLMS','GCIPNLMS','Location','northeast','Orientation','vertical')
 
 
%%
figure
plot(1:N,abs(e_NLMS).^2,'LineWidth',2)
hold on
plot(1:N,abs(e_PNLMS).^2,'LineWidth',2)
plot(1:N,abs(e_IPNLMS).^2,'LineWidth',2)
plot(1:N,abs(e_GCIPNLMS).^2,'LineWidth',2)
 
title('Errors of All Algorithms')
grid on
hold off;
xlabel('n')
ylabel('e[n]^2')
legend('NLMS','PNLMS','IPNLMS','GCIPNLMS','Location','northwest','Orientation','vertical')
 
 


