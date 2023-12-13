%%
clear all;
close all;

%%
% Calculate reference solution for the fast top
! (cd ../.. ; make clean cleantest default test/expandconfig)
! (cd .. ; ./start_test.sh diss_conf.lua)
solcell = load_all_config_and_bin();

%%
% Calculate step sizes and errors
for istep=1:length(solcell)
    solcell{istep}.h = (solcell{istep}.te - solcell{istep}.t0)/solcell{istep}.steps;
end

refpattern.steps = 2*2^21;
solcell = calc_errors(solcell, refpattern);

%%
% Plotting
yyyplot = {'err.abs.q', 'err.abs.v', 'err.abs.l', 'err.abs.lnh'};
xplot = 'h';

makexyyyplot(solcell, [], xplot, yyyplot)
matlab2csv('../../../out/rolling_disk_conv/')

figure;
snapshotplot(solcell{1}, [0]);
title('snapshotplot')
legend('middle trajectory','point trajectory','disk')
view([0, 90]);
matlab2csv('../../../out/rolling_disk_snapshots')
