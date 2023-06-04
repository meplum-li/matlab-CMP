clearvars
tic
mu0 = 0.6;
Ui = 1.74912443;
Sample = parameter();
options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','iter');
options.FunctionTolerance=1e-14;
options.OptimalityTolerance = 1e-8;
options.SpecifyObjectiveGradient = true;
options.SubproblemAlgorithm="factorization";
[delta, fval] = lsqnonlin(@(x) diff_delta_jaco(x),0.15,[],[],options)
toc