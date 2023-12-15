function solcell = create_h(solcell)

[n,m] = size(solcell);

for i = 1:max(n,m)
    if solcell{i}.step_control == 0
        solcell{i}.h = (solcell{i}.te-solcell{i}.t0)/solcell{i}.steps;
    else
        nsteps = length(solcell{i}.rslt.t);
        solcell{i}.h = zeros(nsteps,1);
        for j = 1:nsteps-1
            solcell{i}.h(j) = solcell{i}.rslt.t(j+1)-solcell{i}.rslt.t(j);
        end
        disp('Vector of time step sizes created!')
    end
    
end

end