function Sample = parameter()
Sample.mu=0.6;
Sample.h = 0.;
Sample.delta = 0.3;
Sample.alphaR = 0.;
Sample.periodicity = 1;%1D case
Sample.eta = 1E-6;

%%%for mu
Sample.A_mu_exU = 0.8;%absolute value
Sample.A_mu_exD = 0.6;
Sample.A_mu_exBP = 0.7;
Sample.gammaU = 2*pi*1E-6;%上导线的gamma
Sample.gammaD = 2*pi*1E-4;%上导线的gamma
% gammaSC = 2*pi*1E-1;%下导线的gamma
Sample.gammaBP = 2*pi*1E-4*0;%虚拟导线的gamma
end