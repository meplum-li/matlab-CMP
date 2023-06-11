clearvars
tic
mu0=0.632896719558949;
Sample = parameter();
current_RelTol = Sample.current.RelTol;
current_AbsTol = Sample.current.AbsTol;
% int=integral(@(EF) i_E(Sample,mu0,EF),-inf,inf,"ArrayValued",true,'RelTol',current_RelTol,'AbsTol',current_AbsTol);
ub = max(abs(Sample.A_mu_exU-mu0),abs(Sample.A_mu_exD-mu0)); 
int=integral(@(EF) i_E(Sample,EF, mu0),-ub,+ub,"ArrayValued",true,'RelTol',current_RelTol,'AbsTol',current_AbsTol);
sum(int)
toc
%%
function i_vec = i_E(Sample,EF,mu0)
h = Sample.h;
delta = Sample.delta;
alphaR = Sample.alphaR;
periodcity = Sample.periodicity;
eta = Sample.eta;
sigma0 = eye(2);
sigmaX=[0,1;1,0];
sigmaY=[0,-1i;1i,0];
sigmaZ=[1,0;0,-1];
T_0 = (2-mu0)*kron(sigmaZ, eye(2)) + h*kron(sigmaZ,sigmaZ) - kron( real(delta)*sigmaY+imag(delta)*sigmaX,sigmaY );
T_x = -1*kron(sigmaZ, eye(2)) + alphaR/(2i)*kron(sigmaZ, sigmaY);

N_cen = 30;%中心区长度
%%%电极设置
%都是金属
A_mu_exU = Sample.A_mu_exU;%absolute value
A_mu_exD = Sample.A_mu_exD;
A_mu_exBP = Sample.A_mu_exBP;
mu_exU = A_mu_exU - mu0;%relative value
mu_exD = A_mu_exD - mu0;
mu_exBP = A_mu_exBP - mu0;
gammaU = Sample.gammaU;%上导线的gamma
gammaD = Sample.gammaD;%上导线的gamma
gammaBP = Sample.gammaBP;
%%wide-band limit of real leads
SigmaU = kron(speye(N_cen), -1i*gammaU/2*speye(4));%左导线耦合到中心区的自能
SigmaD = kron(speye(N_cen), -1i*gammaD/2*speye(4));%导线自能
SigmaBP = kron(speye(N_cen), -1i*gammaBP/2*speye(4));
GammaU = 1i*(SigmaU-SigmaU');%左导线的线宽函数
GammaD = 1i*(SigmaD-SigmaD');%右导线的线宽函数
GammaBP = 1i*(SigmaBP - SigmaBP');

H = kron(speye(N_cen), T_0) + kron(diag(ones(N_cen-1,1), 1), T_x) + kron(diag(ones(N_cen-1,1), -1), T_x');
if periodcity == 1
    H = H + kron( diag(1, N_cen-1), T_x ) + kron( diag(1, -N_cen+1), T_x' );
end
GR = inv(eye(4*N_cen)*EF - H-SigmaU - SigmaD - SigmaBP);
SigmaUless = 1i*repmat([EF<mu_exU; EF<mu_exU;EF<-mu_exU; EF<-mu_exU], N_cen, 1).*GammaU;
SigmaDless = 1i*repmat([EF<mu_exD; EF<mu_exD;EF<-mu_exD; EF<-mu_exD], N_cen, 1).*GammaD;
SigmaBPless = 1i*repmat([EF<mu_exBP; EF<mu_exBP;EF<-mu_exBP; EF<-mu_exBP], N_cen, 1).*GammaBP;
Sigmaless= SigmaUless + SigmaDless + SigmaBPless;
Gless = GR * Sigmaless * GR';

%%%计算上下普通电极的电流
i_vec(1) = real(trace( ...
    Gless * ( SigmaU'.*repmat([1,1,-1,-1],1,N_cen)-repmat([1;1;-1;-1],N_cen,1).*SigmaU ) ...
    + GR * (SigmaUless.*repmat([1,1,-1,-1],1,N_cen)) - GR' * (repmat([1;1;-1;-1],N_cen,1).*SigmaUless) ...
    ));
i_vec(2) = real(trace( ...
    Gless * ( SigmaD'.*repmat([1,1,-1,-1],1,N_cen)-repmat([1;1;-1;-1],N_cen,1).*SigmaD ) ...
    + GR * (SigmaDless.*repmat([1,1,-1,-1],1,N_cen)) - GR' * (repmat([1;1;-1;-1],N_cen,1).*SigmaDless) ...
    ));
i_vec(3) = real(trace( ...
    Gless * ( SigmaBP'.*repmat([1,1,-1,-1],1,N_cen)-repmat([1;1;-1;-1],N_cen,1).*SigmaBP ) ...
    + GR * (SigmaBPless.*repmat([1,1,-1,-1],1,N_cen)) - GR' * (repmat([1;1;-1;-1],N_cen,1).*SigmaBPless) ...
    ));
end