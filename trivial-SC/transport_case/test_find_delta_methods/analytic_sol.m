%%% 能隙方程的解析形式，用数值方法在k空间积分，取得吸引势的理论值
clearvars
int=integral(@(k) intg(k),-pi,pi);
Ui = real(1./(int/(4*pi)));

function result = intg(k)
Sample = parameter();
%%%能带结构
mu0=Sample.mu;
h = Sample.h;
delta = Sample.delta;
alphaR = Sample.alphaR;
periodcity = Sample.periodicity;
eta = Sample.eta*0;

result = 1./sqrt( (2-2*cos(k)-mu0-1i*eta).^2 + delta^2 );
end