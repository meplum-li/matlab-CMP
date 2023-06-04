clearvars
tic

mu0 = 0.6;
Ui = 1.74912443;
% Ui = 1.7480;
Sample = parameter();
% delta_curr = delta0*sqrt(1-(Sample.gammaU+Sample.gammaD)/delta0)+0.1;
% delta_next = zero_bias_delta(Ui, mu0, abs(delta_curr));
% options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','iter');
% options.FunctionTolerance=1e-14;
% options.OptimalityTolerance = 1e-8;
% options.SubproblemAlgorithm="factorization";
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
% options.Algorithm = 'levenberg-marquardt';
% [delta,fval]=fsolve(@(delta) zero_bias_delta(delta),0.3)
% [delta, fval] = lsqnonlin(@(delta) zero_bias_delta(delta),0.3,[],[],options)
delta = fmincon(@(delta) zero_bias_delta(delta),0.3,[],[],[],[],[],[],[],options)
% [delta, fval] = fzero(@(delta) real(zero_bias_delta(delta)),0.3)
toc