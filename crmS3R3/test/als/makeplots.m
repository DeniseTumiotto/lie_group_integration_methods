%% %%%%
clear pattern refpattern
pattern.n = 8;

refpattern = pattern;
refpattern.steps = 15 * 2^16;
%%%
solcell = makeplot(pattern, refpattern, 'abs');

%% %%%%
clear pattern refpattern
pattern.N = 2^4;
pattern.banded_iteration_matrix = 1;
pattern.stab2 = 1;

refpattern = pattern;
refpattern.steps = 15 * 2^15;
%%%
solcell = makeplot(pattern, refpattern, 'abs', 'cpu');

%% %%%%
clear pattern refpattern
pattern.N = 2^5;
pattern.banded_iteration_matrix = 0;
pattern.stab2 = 1;

refpattern = pattern;
refpattern.steps = 15 * 2^15;
%%%
solcell = makeplot(pattern, refpattern, 'abs', 'cpu');

%% %%%%
clear pattern refpattern
pattern.N = 2^4;
pattern.stab2 = 0;
pattern.no_K = 0;

refpattern.N = 2^4;
refpattern.stab2 = 1;
refpattern.banded_iteration_matrix = 1;
refpattern.steps = 15 * 2^15;
%%%
solcell = makeplot(pattern, refpattern, 'abs', 'cpu');

% %% %%%%
% clear pattern refpattern
% pattern.stab2 = 0;
% pattern.no_K = 1;
% 
% refpattern.stab2 = 1;
% refpattern.banded_iteration_matrix = 1;
% refpattern.steps = 15 * 2^15;
% %%%
% solcell = makeplot(pattern, refpattern, 'abs', 'cpu');
