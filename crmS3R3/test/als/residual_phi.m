function sols = residual_phi(sols)

n = max(size(sols));

for i = 1:n
    sols{i}.rslt.norm_Phi = vecnorm(sols{i}.rslt.Phi);
    sols{i}.rslt.max_Phi = max(sols{i}.rslt.norm_Phi);
end

end