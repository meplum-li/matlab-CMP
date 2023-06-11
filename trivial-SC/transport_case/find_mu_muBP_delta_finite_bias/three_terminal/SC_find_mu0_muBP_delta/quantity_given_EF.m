function [Tmatrix, GRGA41] = quantity_given_EF(Sample, EF,mu0, A_mu_exBP, delta)
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
N_cen = 30;
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

GR= inv(eye(4*N_cen)*EF - H-SigmaU - SigmaD - SigmaBP);
%%% 计算transmission 矩阵
GAMMA_T = [GammaU; GammaD; GammaBP];%for calculating transmission
temT = GAMMA_T*GR*[GammaU, GammaD, GammaBP].*conj( kron(ones(3),GR) );%稀疏矩阵的写法提速10倍，非稀疏的就不用写成稀疏矩阵了
% %保留自旋自由度
% sumM=blkdiag(eye(8));
% T = real(sumM*temT*sumM');
%自旋缩并
sumM=kron(speye(2*3*N_cen), ones(1,2));
temTT = real(sumM*temT*sumM');
%缩并电极内部指标，只剩下up、down、BP以及各自的电子空穴指标
sumM=kron(  speye(3), kron( ones(1, N_cen), speye(2) )  );
Tmatrix= real(sumM*temTT*sumM');
GRGA41 = GR(4,:)*GR(1,:)';%本身就是实数，加个real去除0i的虚部，便于后续计算
end