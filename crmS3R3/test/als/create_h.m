function solcell = create_h(solcell)

[n,m] = size(solcell);

for i = 1:max(n,m)
    solcell{i}.h = (solcell{i}.te-solcell{i}.t0)/solcell{i}.steps;
end

end