clearvars
tic
Sample = parameter();
options = optimoptions(@fsolve,'Algorithm','trust-region-dogleg','Display','iter');
options.FunctionTolerance=1e-14;
options.OptimalityTolerance = 1e-8;
options.SpecifyObjectiveGradient = true;
options.SubproblemAlgorithm="factorization";
[x, fval] = fsolve(@(x) current_mu_muBP_Jacobi(x),[0.63,0.63],options)
toc