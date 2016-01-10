function [val] = p_vq1q2(v, q1, q2, mu_vq1, mu_vq2,sigma_vq1, sigma_vq2)
    val = normpdf(v, mu_vq1(q1) + mu_vq2(q2), sqrt(sigma_vq1(q1)^2 + sigma_vq2(q2)^2));
end