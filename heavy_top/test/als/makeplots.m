%% %%%%
%problemsets = [30 50 60 70] + 3;
problemsets = 20;

for ii = 1:length(problemsets)
   pattern.problemset = problemsets(ii);
   %pattern.steps = 2^22;

   %refpattern.problemset = problemsets(ii);
   %refpattern.problemset = 33;
   refpattern.steps = 10*2^18;

   %%%
   %solcell = makeplot(pattern, refpattern, 'abs','cpu','SO3xR3');
   solcell = makeplot(pattern, refpattern, 'abs','cpu');
   clear pattern refpattern;
end

swapplots