function mu = bisection(PHImm, Amm, Pk, mu_low, mu_high, tol, max_iter)
iter = 0;
while(1)
    iter = iter + 1;
    mu_mid = (mu_low + mu_high) / 2;
    SUM = 0;
    for m = 1:length(Amm)
        SUM = SUM + (PHImm(m)/(Amm(m) + mu_mid)^2);
    end
    f_mu = SUM - Pk;
    if f_mu > 0
        mu_low = mu_mid;
    else
        mu_high = mu_mid;
    end
    if (abs(mu_high - mu_low) < tol) || (iter == max_iter)
        break;
    end
mu = (mu_low + mu_high)/2;
end