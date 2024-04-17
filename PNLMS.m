function [e_PNLMS,y,h_adapt_final_PNLMS] = PNLMS(d,x,mu_PNLMS,M)
    N = length(x);
    mu_PNLMS = 0.01;
    h_adapt_PNLMS = ones(1,M);
    sigma_PNLMS = 0.01;
    ro_PNLMS = 0.01;
    sigma_V_e = 1;
    V_e = sigma_V_e*randn(N,1);
    Epsilon_PNLMS = 0.00001;
    for n = 1:N-M
        y = h_adapt_PNLMS(n,:)*x(n:n+M-1); % filter coefficients are reversed
        e_PNLMS(n) = d(n+M-1) - y;
        AllFilterCoefficients_PNLMS(n, :) = h_adapt_PNLMS(n,:);
        %Calculation of g and g_bar
        l_inf = max(abs(h_adapt_PNLMS(n,:)));
        l_prime_inf = max(sigma_PNLMS, l_inf);

        for gg = 1:M
            g_PNLMS(gg) = max(ro_PNLMS*l_prime_inf, abs(h_adapt_PNLMS(n,gg)));
        end

        g_bar_PNLMS = (1/M)*sum(g_PNLMS);
        AlgorithmInput(n) = V_e(n) + e_PNLMS(n);
        h_adapt_PNLMS(n+1,:) = h_adapt_PNLMS(n,:) + (flip(g_PNLMS)/(Epsilon_PNLMS +g_bar_PNLMS)).*((mu_PNLMS /(Epsilon_PNLMS + norm(x(n:n+M-1))^2))*AlgorithmInput(n)*x(n:n+M-1)');
    end
        h_adapt_final_PNLMS = h_adapt_PNLMS(end,:);
        h_adapt_final_PNLMS = fliplr(h_adapt_final_PNLMS);
end