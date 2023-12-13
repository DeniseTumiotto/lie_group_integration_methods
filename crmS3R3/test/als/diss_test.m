%%
clear all;
close all;

%%% First, the tests for constant step sizes

%%
% Calculate reference solution
! (cd ../.. ; make clean cleantest default test/expandconfig INTEGRATOR=RATTLie INTEGRATORP=../RATTLie RELAXATION=1 USE_INDEX_1_APPROX=1 TWO_ITERATIONS_MINIMUM=1)
! (cd .. ; ./start_test.sh diss_test_fs_time_ref.lua)
ref = load_latest_config_and_bin();

%%
% Calculate time tests with RATTLie
! (cd ../.. ; make clean cleantest default test/expandconfig INTEGRATOR=RATTLie INTEGRATORP=../RATTLie)
! (cd .. ; ./start_test.sh diss_test_fs_time.lua)
RL = load_all_config_and_bin();

%%
% Calculate time tests with SHAKELie
! (cd ../.. ; make clean cleantest default test/expandconfig INTEGRATOR=SHAKELie INTEGRATORP=../SHAKELie)
! (cd .. ; ./start_test.sh diss_test_fs_time.lua)
SH = load_all_config_and_bin();

%%
% Calculate time tests with gena
! (cd ../.. ; make clean cleantest default test/expandconfig INTEGRATOR=gena INTEGRATORP=../gena)
! (cd .. ; ./start_test.sh diss_test_fs_time.lua)
gena = load_all_config_and_bin();

%%
% Calculate time tests with BLieDF
! (cd ../.. ; make clean cleantest default test/expandconfig INTEGRATOR=BLieDF INTEGRATORP=../BLieDF)
! (cd .. ; ./start_test.sh diss_test_fs_time.lua)
BDF = load_all_config_and_bin();

%%
% Calculate time tests with RATTLie and index-1 approximation
! (cd ../.. ; make clean cleantest default test/expandconfig INTEGRATOR=RATTLie INTEGRATORP=../RATTLie USE_INDEX_1_APPROX=1)
! (cd .. ; ./start_test.sh diss_test_fs_time.lua)
RL_i1 = load_all_config_and_bin();

%%
% Calculate time tests with RATTLie with only lm and lp
! (cd ../.. ; make clean cleantest default test/expandconfig INTEGRATOR=RATTLie INTEGRATORP=../RATTLie ONLY_LM_LP=1)
! (cd .. ; ./start_test.sh diss_test_fs_time.lua)
RL_lmp = load_all_config_and_bin();

%%
save diss_test_time_novar;

%%
% Calculate errors
ref.external =  [ref.external '_ref'];
refpattern.external = ref.external;
refpattern.steps = ref.steps;
refpattern.integrator = ref.integrator;

RL = calc_errors({RL{:}, ref}, refpattern);
SH = calc_errors({SH{:}, ref}, refpattern);
gena = calc_errors({gena{:}, ref}, refpattern);
BDF = calc_errors({BDF{:}, ref}, refpattern);
RL_i1 = calc_errors({RL_i1{:}, ref}, refpattern);
RL_lmp = calc_errors({RL_lmp{:}, ref}, refpattern);

%%
% Add step sizes
for i=1:length(RL)
   RL{i}.h = (RL{i}.te - RL{i}.t0)/RL{i}.steps;
end
for i=1:length(SH)
   SH{i}.h = (SH{i}.te - SH{i}.t0)/SH{i}.steps;
end
for i=1:length(gena)
   gena{i}.h = (gena{i}.te - gena{i}.t0)/gena{i}.steps;
end
for i=1:length(BDF)
   BDF{i}.h = (BDF{i}.te - BDF{i}.t0)/BDF{i}.steps;
end
for i=1:length(RL_i1)
   RL_i1{i}.h = (RL_i1{i}.te - RL_i1{i}.t0)/RL_i1{i}.steps;
end
for i=1:length(RL_lmp)
   RL_lmp{i}.h = (RL_lmp{i}.te - RL_lmp{i}.t0)/RL_lmp{i}.steps;
end


%%
% Make and save plots (don't plot the reference solution)
makexyyyplot(RL(1:end-1), [], 'h', {'err.abs.q','err.abs.v','err.abs.l'});
matlab2csv('../../../out/crm_fs_RATTLie_abserr/');

makexyyyplot(SH(1:end-1), [], 'h', {'err.abs.q','err.abs.v','err.abs.l'});
matlab2csv('../../../out/crm_fs_SHAKELie_abserr/');

makexyyyplot(gena(1:end-1), [], 'h', {'err.abs.q','err.abs.v','err.abs.l'});
matlab2csv('../../../out/crm_fs_gena_abserr/');

makexyyyplot(BDF(1:end-1), [], 'h', {'err.abs.q','err.abs.v','err.abs.l'});
matlab2csv('../../../out/crm_fs_BLieDF_abserr/');

makexyyyplot(RL_i1(1:end-1), [], 'h', {'err.abs.q','err.abs.v','err.abs.l'});
matlab2csv('../../../out/crm_fs_RATTLie_i1_abserr/');

makexyyyplot(RL_lmp(1:end-1), [], 'h', {'err.abs.q','err.abs.v','err.abs.l','err.abs.lm','err.abs.lmp'});
matlab2csv('../../../out/crm_fs_RATTLie_lmp_abserr/');

%%
%%% Variable step sizes

%% Calculate reference solution for the variable steps size scenario
! (cd ../.. ; make clean cleantest default test/expandconfig INTEGRATOR=RATTLie INTEGRATORP=../RATTLie RELAXATION=1 USE_INDEX_1_APPROX=1 TWO_ITERATIONS_MINIMUM=1)
! (cd .. ; ./start_test.sh diss_test_fs_time_var_ref.lua)
ref_var = load_latest_config_and_bin();

%%
% Calculate time tests with RATTLie with variable steps
! (cd ../.. ; make clean cleantest default test/expandconfig INTEGRATOR=RATTLie INTEGRATORP=../RATTLie VARIABLE_STEPS=1)
! (cd .. ; ./start_test.sh diss_test_fs_time_var.lua)
RL_var = load_all_config_and_bin();


%%
% Calculate time tests with RATTLie with variable steps and index-1 approximation
! (cd ../.. ; make clean cleantest default test/expandconfig INTEGRATOR=RATTLie INTEGRATORP=../RATTLie VARIABLE_STEPS=1 USE_INDEX_1_APPROX=1)
! (cd .. ; ./start_test.sh diss_test_fs_time_var.lua)
RL_var = load_all_config_and_bin();

%%
save diss_test_time_var;

%%
% Calculate errors
ref.external =  [ref.external '_ref'];
refpattern.external = ref.external;
refpattern.steps = ref.steps;
refpattern.integrator = ref.integrator;

RL_var = calc_interp_errors(RL_var, {ref_var}, []);
RL_var_i1 = calc_interp_errors(RL_var_i1, {ref_var}, []);

%%
% Add step sizes
for i=1:length(RL_var)
   RL_var{i}.hmax = max(RL_var{i}.rslt.t(2:end) - RL_var{i}.rslt.t(1:end-1));
end
for i=1:length(RL_var_i1)
   RL_var_i1{i}.hmax = max(RL_var_i1{i}.rslt.t(2:end) - RL_var_i1{i}.rslt.t(1:end-1));
end


%%
% Make and save plots
makexyyyplot(RL_var, [], 'hmax', {'err.abs.q','err.abs.v','err.abs.l'});
matlab2csv('../../../out/crm_fs_RATTLie_var_abserr/');

makexyyyplot(RL_var_i1, [], 'hmax', {'err.abs.q','err.abs.v','err.abs.l'});
matlab2csv('../../../out/crm_fs_RATTLie_var_i1_abserr/');
