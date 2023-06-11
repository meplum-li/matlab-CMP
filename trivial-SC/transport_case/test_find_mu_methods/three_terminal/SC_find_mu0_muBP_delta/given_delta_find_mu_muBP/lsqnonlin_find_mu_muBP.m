clearvars
tic
Sample = parameter();
options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','iter');
options.FunctionTolerance=1e-14;
options.OptimalityTolerance = 1e-8;
% options.StepTolerance=1e-9;
options.SpecifyObjectiveGradient = true;
options.SubproblemAlgorithm="factorization";
[mu, fval] = lsqnonlin(@(x) current_mu_muBP_Jacobi(x),[0.63,0.63],[],[],options)
toc