clearvars
tic
mu0 = 0.6;
Ui = 1.74912443;
Sample = parameter();
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
options.Algorithm = "trust-region-reflective";
options.SpecifyObjectiveGradient =true;
options.HessianFcn = 'objective';
delta = fmincon(@(x) diff_deltaSQ_jaco_Hess(x),0.15,[],[],[],[],[],[],[],options)
toc