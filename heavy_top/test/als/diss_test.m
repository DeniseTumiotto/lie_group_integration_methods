%%
clear all;
close all;

%%
% Calculate reference solution for the fast top
! (cd ../.. ; make clean cleantest gena default test/expandconfig)
! (cd .. ; ./start_test.sh diss_config_ref_fast.lua)
reffast = load_latest_config_and_bin();

% fast top: do it one more time!
for fastit=1:100

%%
% Calculate tests for the fast top
! (cd ../.. ; make clean cleantest RATTLie default test/expandconfig)
! (cd .. ; ./start_test.sh diss_config_test_fast.lua)

%%
% Load results for the fast top
liegroups = 3:7;
for ilie=1:length(liegroups)
    pattern.liegroup = liegroups(ilie);
    fast{ilie} = load_all_config_and_bin(pattern);
    [~,ind] = sort(catsolcell(fast{ilie},'steps'));
    fast{ilie}(:) = fast{ilie}(ind);
end

%%
% Calculate step sizes, errors and relative computation times
reflieind = find(liegroups == 5);
for ilie=1:length(liegroups)
    for istep=1:length(fast{ilie})
        fast{ilie}{istep}.h = (fast{ilie}{istep}.te - fast{ilie}{istep}.t0)/fast{ilie}{istep}.steps;

        ind = find(catsolcell(fast{reflieind},'steps') == fast{ilie}{istep}.steps);
        if length(ind) ~= 1
            fast{ilie}{istep}.rel_cpu = NaN;
        else
            fast{ilie}{istep}.rel_cpu = fast{ilie}{istep}.stats.cpu_time/fast{reflieind}{ind}.stats.cpu_time;
        end
    end

    fast{ilie} = calc_SO3xR3_errors(fast{ilie}, {reffast})

    fastrelcpu(fastit,ilie,:) = catsolcell(fast{ilie},'rel_cpu');
end
end %for fastit=1:?

%%
% add the mean of the fastrelcpu to te fast solcell
for ilie=1:length(liegroups)
    for istep=1:length(fast{ilie})
        fast{ilie}{istep}.avg_rel_cpu = mean(fastrelcpu(:,ilie,istep));
    end
end

%%
% Calculate reference solution for the slow top
! (cd ../.. ; make clean cleantest gena default test/expandconfig)
! (cd .. ; ./start_test.sh diss_config_ref_slow.lua)
refslow = load_latest_config_and_bin();

%%
% Calculate tests for the slow top
! (cd ../.. ; make clean cleantest RATTLie default test/expandconfig)
! (cd .. ; ./start_test.sh diss_config_test_slow.lua)


%%
% Load results for the slow top
liegroups = 3:7;
for ilie=1:length(liegroups)
    pattern.liegroup = liegroups(ilie);
    slow{ilie} = load_all_config_and_bin(pattern);
    [~,ind] = sort(catsolcell(slow{ilie},'steps'));
    slow{ilie}(:) = slow{ilie}(ind);
end

%%
% Calculate step sizes and relative computation times
reflieind = find(liegroups == 5);
for ilie=1:length(liegroups)
    for istep=1:length(slow{ilie})
        slow{ilie}{istep}.h = (slow{ilie}{istep}.te - slow{ilie}{istep}.t0)/slow{ilie}{istep}.steps;
    end

    slow{ilie} = calc_SO3xR3_errors(slow{ilie}, {refslow});
end


save('disstest.mat')

%%
%'SO3xR3_err.abs.q'
%makexyplot(sc, [], 'steps', 'SO3xR3_err.abs.q','liegroup')
liegroupstoplot = find(ismember(liegroups,liegroups));
makexyplot([fast{liegroupstoplot}], [], 'h', 'avg_rel_cpu','liegroup')
matlab2csv('../../../out/heavy_top_fast_rel_cpu/')
close gcf;

liegroupstoplot = find(ismember(liegroups,liegroups));
makexyplot([fast{liegroupstoplot}], [], 'h', 'SO3xR3_err.abs.q','liegroup');
matlab2csv('../../../out/heavy_top_fast_abs_err/')
close gcf;

liegroupstoplot = find(ismember(liegroups,liegroups));
makexyplot([slow{liegroupstoplot}], [], 'h', 'SO3xR3_err.abs.q','liegroup');
matlab2csv('../../../out/heavy_top_slow_abs_err/')
close gcf;

%%
close all;
clear all;

%%
% Make energy plots
! (cd ../.. ; make clean cleantest gena default test/expandconfig)
! (cd .. ; ./start_test.sh diss_config_test_energy.lua)
engena = load_latest_config_and_bin();

! (cd ../.. ; make clean cleantest RATTLie default test/expandconfig)
! (cd .. ; ./start_test.sh diss_config_test_energy.lua)
enRL = load_latest_config_and_bin();

! (cd ../.. ; make clean cleantest BLieDF default test/expandconfig)
! (cd .. ; ./start_test.sh diss_config_test_energy.lua)
enBDF = load_latest_config_and_bin();

save('disstest_energy.mat')

%%
% Calculate energies
engena = calc_energy(engena);
enRL = calc_energy(enRL);
enBDF = calc_energy(enBDF);


%%
% actual plot
figure();
plot(engena.rslt.t, engena.rslt.energy,...
     enRL.rslt.t, enRL.rslt.energy,...
     enBDF.rslt.t, enBDF.rslt.energy);
legend('gena','RATTLie','BLieDF')
title('mechanic energy over time')
xlabel('t')
ylabel('energy')
matlab2csv('../../../out/heavy_top_fast_energy/');
