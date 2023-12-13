function sol = calc_energy(sol)

for i=1:length(sol.rslt.t)
   sol.rslt.energy_kinetic(i) = 1/2*sol.rslt.v(:,i)'*([sol.inerJ,sol.mass,sol.mass,sol.mass]'.*sol.rslt.v(:,i));
   sol.rslt.energy_potential(i) = -sol.mass * sol.gravity*sol.rslt.q(5:7,i);
end

sol.rslt.energy = sol.rslt.energy_kinetic + sol.rslt.energy_potential;