clearvars
tic
mu0 = 0.6;
% Ui = 1.74912443;
Sample = parameter();
% Sample.gammaU = 0.152;%上导线的gamma
% Sample.gammaD = 0.152;%上导线的gamma
options = optimoptions(@fsolve,'Algorithm','trust-region-dogleg','Display','iter');
options.FunctionTolerance=1e-14;
options.OptimalityTolerance = 1e-10;
options.SpecifyObjectiveGradient = true;
options.SubproblemAlgorithm="factorization";
options.UseParallel = true;
[delta, fval] = fsolve(@(x) diff_delta_jaco(Sample, mu0, x),2,options)
toc