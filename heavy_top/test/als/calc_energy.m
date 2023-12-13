function sol = calc_energy(sol)
% Calculate energies

if not(sol.liegroup == 6)
   error('Lie group needs to be 6 (S^3 |x R^3)');
end

if not(sol.rslt.finished == 1)
   warning('solution did not finish');
end

sol.rslt.energy_kinetic = zeros(1,length(sol.rslt.t))
sol.rslt.energy_potential = zeros(1,length(sol.rslt.t))

for it = 1:length(sol.rslt.t)
   sol.rslt.energy_kinetic(it) = (sol.rslt.v(1:3,it).*sol.inerJ')' * sol.rslt.v(1:3,it)/2 ...
                               +  sol.mass*sol.rslt.v(4:6,it)'*sol.rslt.v(4:6,it)/2;
   sol.rslt.energy_potential(it) = -sol.mass*sol.gravity*sol.rslt.q(5:7,it);
end

sol.rslt.energy = sol.rslt.energy_kinetic + sol.rslt.energy_potential;
sol.rslt.energy = sol.rslt.energy - sol.rslt.energy(1);
