clearvars
tic
mu0 = 0.6;
Ui = 1.74912443;
Sample = parameter();
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
% options.Algorithm = "trust-region-reflective";
% options.SpecifyObjectiveGradient =true;
% delta = fmincon(@(x) diff_delta_jaco_Hess(x),0.15,[],[],[],[],[],[],[],options)
delta = fmincon(@(x) zero_bias_delta(x),0.15,[],[],[],[],[],[],[],options)

toc