clearvars
tic
Sample = parameter();
options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','iter');
options.FunctionTolerance=1e-10;
options.OptimalityTolerance = 1e-6;
% options.StepTolerance=1e-9;
options.SpecifyObjectiveGradient = true;
options.SubproblemAlgorithm="factorization";
[x, fval] = lsqnonlin(@(x) mu_muBP_Delta_Jacobi(Sample, x),[0.63,0.63,2],[],[],options);
% [x, fval] = lsqnonlin(@(x) v2_mu_muBP_Delta_Jacobi(Sample, x),[0.63,0.63,2],[],[],options);
% [x, fval] = lsqnonlin(@(x) v3_mu_muBP_Delta_Jacobi(Sample, x),[0.63,0.63,2],[],[],options);
% [x, fval] = lsqnonlin(@(x) v4_mu_muBP_Delta_Jacobi(Sample, x),[0.65,0.65,2],[],[],options);
% [x, fval] = lsqnonlin(@(x) v5_mu_muBP_Delta_Jacobi(Sample, x),[0.65,0.65,2],[],[],options);
toc