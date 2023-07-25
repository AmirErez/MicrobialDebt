function S=calc_entropy_nats(rho)
% Returns the entropy in nats


S = NaN*zeros(size(rho, 2),1);
for ii=1:size(rho, 2)
    norm_b = rho(:,ii)/sum(rho(:,ii));
    lognorm = log(norm_b);
    lognorm(lognorm==-Inf) = 0;
    S(ii) = -sum(norm_b.*lognorm);
end
