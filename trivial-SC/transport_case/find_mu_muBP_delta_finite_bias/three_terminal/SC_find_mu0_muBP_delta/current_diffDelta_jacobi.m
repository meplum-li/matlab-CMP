function integrand = current_diffDelta_jacobi(Sample, EF, mu0, A_mu_exBP, delta)
h = Sample.h;
alphaR = Sample.alphaR;
periodcity = Sample.periodicity;
sigma0 = eye(2);
sigmaX=[0,1;1,0];
sigmaY=[0,-1i;1i,0];
sigmaZ=[1,0;0,-1];
T_0 = (2-mu0)*kron(sigmaZ, eye(2)) + h*kron(sigmaZ,sigmaZ) - kron( real(delta)*sigmaY+imag(delta)*sigmaX,sigmaY );
T_x = -1*kron(sigmaZ, eye(2)) + alphaR/(2i)*kron(sigmaZ, sigmaY);

%% transport
N_cen = Sample.N_cen;
%%%电极设置
A_mu_exU = Sample.A_mu_exU;%absolute value
A_mu_exD = Sample.A_mu_exD;
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
H_partial_mu =kron( speye(N_cen), -kron(sigmaZ, eye(2)) );
H_partial_delta =kron(speye(N_cen), -kron(sigmaY, sigmaY ));
% GR= (eye(4*N_cen)*EF - H-SigmaU - SigmaD - SigmaBP)\eye(4*N_cen);
GR= inv(eye(4*N_cen)*EF - H-SigmaU - SigmaD - SigmaBP);
GR_partial_mu = GR * H_partial_mu * GR;
GR_partial_delta = GR*H_partial_delta*GR;
SigmaUless = 1i*repmat([EF<mu_exU; EF<mu_exU;EF<-mu_exU; EF<-mu_exU], N_cen, 1).*GammaU;
SigmaDless = 1i*repmat([EF<mu_exD; EF<mu_exD;EF<-mu_exD; EF<-mu_exD], N_cen, 1).*GammaD;
SigmaBPless = 1i*repmat([EF<mu_exBP; EF<mu_exBP;EF<-mu_exBP; EF<-mu_exBP], N_cen, 1).*GammaBP;
Sigmaless= SigmaUless + SigmaDless + SigmaBPless;
Gless = GR * Sigmaless * GR';

%%% 计算transmission 矩阵
GAMMA_T = [GammaU; GammaD; GammaBP];%for calculating transmission
temT = GAMMA_T*GR*[GammaU, GammaD, GammaBP].*conj( kron(ones(3),GR) );%稀疏矩阵的写法提速10倍，非稀疏的就不用写成稀疏矩阵了
temT_partial_mu = GAMMA_T*GR_partial_mu*[GammaU, GammaD, GammaBP].*conj( kron(ones(3),GR) );%稀疏矩阵的写法提速10倍，非稀疏的就不用写成稀疏矩阵了
temT_partial_delta = GAMMA_T*GR_partial_delta*[GammaU, GammaD, GammaBP].*conj( kron(ones(3),GR) );%稀疏矩阵的写法提速10倍，非稀疏的就不用写成稀疏矩阵了
%         temT_partial_mu = temT_partial_mu+temT_partial_mu';
% %保留自旋自由度
% sumM=blkdiag(eye(8));
% T = real(sumM*temT*sumM');
%自旋缩并
sumM=kron(speye(2*3*N_cen), ones(1,2));
temTT = real(sumM*temT*sumM');
temTT_partial_mu = real(sumM*temT_partial_mu*sumM');
temTT_partial_delta = real(sumM*temT_partial_delta*sumM');
%缩并电极内部指标，只剩下up、down、BP以及各自的电子空穴指标
sumM=kron(  speye(3), kron( ones(1, N_cen), speye(2) )  );
TT= real(sumM*temTT*sumM');
TT_partial_mu = 2*real(sumM*temTT_partial_mu*sumM');
TT_partial_delta = 2*real(sumM*temTT_partial_delta*sumM');
%%%计算电流
fUe = (EF<mu_exU);
fUh = (EF<-mu_exU);
fDe = (EF<mu_exD);
fDh = (EF<-mu_exD);
fBPe = (EF<mu_exBP);
fBPh = (EF<-mu_exBP);
current = kron(eye(3), [1,1]) * sum( kron(ones(3,1), [1;-1]).*TT.*( repelem([fUe,fUh, fDe, fDh, fBPe, fBPh].', 1, 6)-repelem([fUe,fUh, fDe, fDh, fBPe, fBPh], 6, 1) ),2 );
% diffDelta = diag( kron(eye(N_cen),[0,0,0,1])*Gless*kron(eye(N_cen),[1;0;0;0]) );
diffDelta = real( 1i*Gless(4,1) );%因为周期性边界条件，不同位置处应该是相同的。
current_Dmu = kron(eye(3), [1,1]) * sum( kron(ones(3,1), [1;-1]).*TT_partial_mu.*( repelem([fUe,fUh, fDe, fDh, fBPe, fBPh].', 1, 6)-repelem([fUe,fUh, fDe, fDh, fBPe, fBPh], 6, 1) ),2 );
current_Ddelta = kron(eye(3), [1,1]) * sum( kron(ones(3,1), [1;-1]).*TT_partial_delta.*( repelem([fUe,fUh, fDe, fDh, fBPe, fBPh].', 1, 6)-repelem([fUe,fUh, fDe, fDh, fBPe, fBPh], 6, 1) ),2 );
%%%diffDelta_Dmu
K_mu = GR_partial_mu*Sigmaless*GR';%可以简化
diffDelta_Dmu = real( 1i*( K_mu(4,1)-conj(K_mu(1,4)) ) );
%%%diffDelta_Ddelta
K_delta = GR_partial_delta * Sigmaless * GR';
diffDelta_Ddelta = real( 1i*( K_delta(4,1) - conj(K_delta(1,4)) ) );
integrand = [current;diffDelta;current_Dmu;current_Ddelta;diffDelta_Dmu;diffDelta_Ddelta];
end