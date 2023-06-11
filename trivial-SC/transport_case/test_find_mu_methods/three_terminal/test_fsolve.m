clearvars
tic
mu0 = 0.6;
Sample = parameter();
options = optimoptions(@fsolve,'Algorithm','trust-region-dogleg','Display','iter');
options.FunctionTolerance=1e-14;
options.OptimalityTolerance = 1e-8;
options.SpecifyObjectiveGradient = true;
options.SubproblemAlgorithm="factorization";
% options.UseParallel = true;
[delta, fval] = fsolve(@(x) current_Jaco(x),0.6,options)
toc