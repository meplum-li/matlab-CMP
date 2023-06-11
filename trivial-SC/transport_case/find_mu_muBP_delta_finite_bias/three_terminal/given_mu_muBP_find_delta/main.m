clearvars
% gamma = linspace(0.01,0.16,20);
gamma = 0.15;
N_gamma = length(gamma);
mu0 = 19/30;
Ui = 1.74912443;
delta = zeros(N_gamma,1);

for ii = 1 : N_gamma
    Sample = parameter();
    Sample.gammaD = gamma(ii);
    Sample.gammaU = gamma(ii);
    options = optimoptions(@fsolve,'Algorithm','trust-region-dogleg','Display','iter');
options.FunctionTolerance=1e-14;
options.OptimalityTolerance = 1e-10;
options.SpecifyObjectiveGradient = true;
options.SubproblemAlgorithm="factorization";
% options.UseParallel = true;
delta(ii) = fsolve(@(x) diff_delta_jaco(Sample,Ui, mu0, x),0.15,options);

end