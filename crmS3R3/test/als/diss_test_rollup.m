%%
clear all;
close all;

%%
% Calculate solutions
! (cd ../.. ; make clean cleantest default test/expandconfig INTEGRATOR=gena INTEGRATORP=../gena)
! (cd .. ; ./start_test.sh diss_test_rollup.lua)
sc = load_all_config_and_bin();

save('diss_test_rollup.mat')

%%
% Add step sizes
for i=1:length(sc)
   sc{i}.h = (sc{i}.te - sc{i}.t0)/sc{i}.steps;
end

%%
% Make and save plots
snapshotplot(sc{2}, [0, 0.5, 1, 1.5, 2, 3, 30], 8);
title('snapshotsplot');
matlab2csv('../../../out/crm_rollup_snapshots/');
close gcf;


figure();
leg = {};
for i=1:length(sc)
   sol = sc{i};
   semilogy(sol.rslt.t, sqrt(sum(sol.rslt.q(end-2:end,:).^2,1)) + sqrt(sum(bsxfun(@minus, sol.rslt.q(end-6:end-3,:), [-1/sqrt(2);0;-1/sqrt(2);0]).^2,1)));
   hold on;
   leg{end+1} = ['h=' num2str(sol.h,'%5.2e')];
end
hold off;
ylabel('distance')
xlabel('t')
title('distance over time')
legend(leg{:});
matlab2csv('../../../out/crm_rollup_end_distance/');


figure();
leg = {};
for i=1:length(sc)
   sol = sc{i};
   semilogy(sol.rslt.t, sqrt(sum(sol.rslt.v(end-5:end,:).^2,1)));
   hold on;
   leg{end+1} = ['h=' num2str(sol.h,'%5.2e')];
end
hold off;
ylabel('norm velocity')
xlabel('t')
title('velocity over time')
legend(leg{:});
matlab2csv('../../../out/crm_rollup_end_velocity/');
