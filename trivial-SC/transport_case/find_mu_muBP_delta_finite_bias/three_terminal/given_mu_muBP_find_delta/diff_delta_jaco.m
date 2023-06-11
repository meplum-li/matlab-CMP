function [diff_delta, J] = diff_delta_jaco(Sample, mu0,delta0)
% function [diff_delta] = zero_bias_delta(Ui, mu0, delta0)
%zero_bias_Ui 连接导线、零偏压的情况、给定孤立中心的吸引势，单次迭代计算新的delta
% 上下电极接上中心区，平衡的情况下，考虑能隙方程，来确定吸引势的大小
% 采取周期性边界条件,只需要输入标量Ui
% 确定的吸引势用来在后续讨论中自洽计算gap
% Ui: 在位型有效吸引势
% mu0: 超导的化学势
% delta0: 初始化的试探delta0
%delta_new_i: 周期性边界条件下返回一个标量
% tic
% delta0 = delta00+1i*0;
% Ui=1.74912443;
% mu0 = 0.6;
% Sample = parameter();
Ui = Sample.Ui;
gap_RelTol = Sample.gap.RelTol;
gap_AbsTol = Sample.gap.AbsTol;
int=integral(@(EF) Gless21(Sample, mu0, delta0, EF),-inf,inf,"ArrayValued",true,'RelTol',gap_RelTol,'AbsTol',gap_AbsTol);
diff_delta = real(Ui/(2*pi)*int(1) - delta0);
J = real(Ui/(2*pi)*int(3) - 1);
%%
    function [integrand] = Gless21(Sample, mu0, delta, EF)
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

        % diffDelta = diag( kron(eye(N_cen),[0,0,0,1])*Gless*kron(eye(N_cen),[1;0;0;0]) );
        diffDelta = 1i*Gless(4,1);%因为周期性边界条件，不同位置处应该是相同的。
        %%%diffDelta_Dmu
        K_mu = GR_partial_mu*Sigmaless*GR';%可以简化
        diffDelta_Dmu =  1i*( K_mu(4,1)-conj(K_mu(1,4)) ) ;
        %%%diffDelta_Ddelta
        K_delta = GR_partial_delta * Sigmaless * GR';
        diffDelta_Ddelta = 1i*( K_delta(4,1) - conj(K_delta(1,4)) ) ;
        integrand = [diffDelta, diffDelta_Dmu,diffDelta_Ddelta ];
end
end