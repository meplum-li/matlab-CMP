%%% 孤立的中心区，平衡的情况下，考虑能隙方程，来确定吸引势的大小
clearvars
tic 
Sample = parameter();
delta = Sample.delta;
% eta = Sample.eta;
int=integral(@(EF) Gless21(EF),-4,0.,"ArrayValued",true,'RelTol',1e-4,'AbsTol',1e-13);
Ui = delta./(-int/(2*pi*1i));
toc
% tic
% EF = linspace(-4,4,200);
% GPN = zeros(30,length(EF));
% for ii = 1 : length(EF)
%     GPN(:,ii) = Gless21(EF(ii));
% end
% %%
% figure
% plot(EF,abs(GPN(3,:)), 'k',LineWidth=2)
% % ylim([0,1E-2])
% xlim([EF(1),EF(end)])
% toc
%%
function result = Gless21(EF)
Sample = parameter();
%%%能带结构
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
T_0 = (2-mu0)*kron(sigmaZ, eye(2)) + h*kron(sigmaZ,sigmaZ) - delta* kron(sigmaY, sigmaY);
T_x = -1*kron(sigmaZ, eye(2)) + alphaR/(2i)*kron(sigmaZ, sigmaY);

N_cen = 30;%中心区长度
H = kron(speye(N_cen), T_0) + kron(diag(ones(N_cen-1,1), 1), T_x) + kron(diag(ones(N_cen-1,1), -1), T_x');
if periodcity == 1
    H = H + kron( diag(1, N_cen-1), T_x ) + kron( diag(1, -N_cen+1), T_x' );
end
GR = inv(eye(4*N_cen)*(EF+1i*eta) - H);
Gless = -(EF<=0)*(GR - GR');
% result=diag(kron(eye(N_cen),[0,0,0,1])*Gless);
result=diag( kron(eye(N_cen),[0,0,0,1])*Gless*kron(eye(N_cen),[1;0;0;0]) );
end