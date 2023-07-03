clearvars
tic
Sample = parameter();
options = optimoptions(@fsolve,'Algorithm','trust-region-dogleg','Display','iter');
options.FunctionTolerance=1e-14;
options.OptimalityTolerance = 1e-6;
options.SpecifyObjectiveGradient = true;
options.SubproblemAlgorithm="factorization";
[x, fval] = fsolve(@(x) v5_mu_muBP_Delta_Jacobi(Sample,x),[0.63,0.63,2],options)
toc