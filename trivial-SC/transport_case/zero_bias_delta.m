function [delta_new_i] = zero_bias_delta(Ui, mu0, delta0)
%zero_bias_Ui 连接导线、零偏压的情况、给定孤立中心的吸引势，单次迭代计算新的delta
% 上下电极接上中心区，平衡的情况下，考虑能隙方程，来确定吸引势的大小
% 采取周期性边界条件,只需要输入标量Ui
% 确定的吸引势用来在后续讨论中自洽计算gap
% Ui: 在位型有效吸引势
% mu0: 超导的化学势
% delta0: 初始化的试探delta0
%delta_new_i: 周期性边界条件下返回一个标量
% tic
Sample = parameter();
%     delta = Sample.delta;
% eta = Sample.eta;
gap_RelTol = Sample.gap.RelTol;
gap_AbsTol = Sample.gap.AbsTol;
int=integral(@(EF) Gless21(Sample, mu0, delta0, EF),-inf,inf,"ArrayValued",true,'RelTol',gap_RelTol,'AbsTol',gap_AbsTol);
% Ui = -Sample.delta*2*pi*1i./int;
delta_new_i = mean(-Ui *int/(2*pi*1i));
% ElapsedTime = toc;
% fprintf('%%------------------------\n')
% fprintf('  delta_new_i\t   Time Cost\n')
% fprintf('%12.8f     %8.2f\n',delta_new_i,ElapsedTime)

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
    function result = Gless21(Sample, mu0, delta, EF)
        %     Sample = parameter();
        %     mu0=Sample.mu;
        h = Sample.h;
%         delta = Sample.delta;
        alphaR = Sample.alphaR;
        periodcity = Sample.periodicity;
        eta = Sample.eta;
        sigma0 = eye(2);
        sigmaX=[0,1;1,0];
        sigmaY=[0,-1i;1i,0];
        sigmaZ=[1,0;0,-1];
        T_0 = (2-mu0)*kron(sigmaZ, eye(2)) + h*kron(sigmaZ,sigmaZ) - kron( real(delta)*sigmaY+imag(delta)*sigmaX,sigmaY );
        T_x = -1*kron(sigmaZ, eye(2)) + alphaR/(2i)*kron(sigmaZ, sigmaY);

        N_cen = Sample.N_cen;%中心区长度
        %%%电极设置
        %都是金属
        A_mu_exU = mu0;%absolute value
        A_mu_exD = mu0;
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
        % result=diag(kron(eye(N_cen),[0,0,0,1])*Gless);
        result=diag( kron(eye(N_cen),[0,0,0,1])*Gless*kron(eye(N_cen),[1;0;0;0]) );
    end
end