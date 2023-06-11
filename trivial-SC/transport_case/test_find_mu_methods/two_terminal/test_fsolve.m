clearvars
tic
mu0 = 0.6;
Sample = parameter();
options = optimoptions(@fsolve,'Algorithm','trust-region-dogleg','Display','iter');
% options.FunctionTolerance=1e-14;
options.OptimalityTolerance = 1e-10;
options.SpecifyObjectiveGradient = true;
options.SubproblemAlgorithm="factorization";
options.UseParallel = true;
[delta, fval] = fsolve(@(x) transmission_to_current(x),0.6,options)
toc