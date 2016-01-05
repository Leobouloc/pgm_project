function val = P_vq(v, q, mu_vq, sigma_vq)
    val = normpdf(v, mu_vq(q), sigma_vq(q));
end