clearvars
tic
Sample = parameter();
current_RelTol = Sample.current.RelTol;
current_AbsTol = Sample.current.AbsTol;
int=integral(@(EF) i_E(EF),-inf,inf,"ArrayValued",true,'RelTol',current_RelTol,'AbsTol',current_AbsTol);
sum(int)
toc

function i_vec = i_E(EF)
Sample = parameter();
mu0=Sample.mu;
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
mu_exU = A_mu_exU - mu0;%relative value
mu_exD = A_mu_exD - mu0;
gammaU = Sample.gammaU;%上导线的gamma
gammaD = Sample.gammaD;%上导线的gamma
%wide-band limit of real leads
SigmaU = kron(speye(N_cen), -1i*gammaU/2*speye(4));%左导线耦合到中心区的自能
SigmaD = kron(speye(N_cen), -1i*gammaD/2*speye(4));%导线自能
GammaU= 1i*(SigmaU-SigmaU');%左导线的线宽函数
GammaD = 1i*(SigmaD-SigmaD');%右导线的线宽函数

H = kron(speye(N_cen), T_0) + kron(diag(ones(N_cen-1,1), 1), T_x) + kron(diag(ones(N_cen-1,1), -1), T_x');
if periodcity == 1
    H = H + kron( diag(1, N_cen-1), T_x ) + kron( diag(1, -N_cen+1), T_x' );
end
GR = inv(eye(4*N_cen)*(EF) - H - SigmaU - SigmaD);
SigmaUless = 1i*repmat([EF<=mu_exU; EF<=mu_exU;EF<=-mu_exU; EF<=-mu_exU], N_cen, 1).*GammaU;
SigmaDless = 1i*repmat([EF<=mu_exD; EF<=mu_exD;EF<=-mu_exD; EF<=-mu_exD], N_cen, 1).*GammaD;
Sigmaless= SigmaUless + SigmaDless;
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
end