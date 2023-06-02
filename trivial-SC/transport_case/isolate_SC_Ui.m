function [Ui] = isolate_SC_Ui(mu0)
%isolate_SC_Ui 零偏压的情况下计算
%%% 上下电极接上中心区，平衡的情况下，考虑能隙方程，来确定吸引势的大小
%%% 采取周期性边界条件更接近理论值
%%% 确定的吸引势用来在后续讨论中自洽计算gap
tic
Sample = parameter();
%     delta = Sample.delta;
% eta = Sample.eta;
gap_RelTol = Sample.gap.RelTol;
gap_AbsTol = Sample.gap.AbsTol;
int=integral(@(EF) Gless21(Sample, mu0, EF),-inf,inf,"ArrayValued",true,'RelTol',gap_RelTol,'AbsTol',gap_AbsTol);
Ui = -Sample.delta*2*pi*1i./int;
ElapsedTime = toc;
fprintf('%%------------------------\n')
fprintf('  Ui\t   Time Cost\n')
fprintf('%12.8f     %8.2f\n',mean(Ui),ElapsedTime)
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
    function result = Gless21(Sample, mu0, EF)
        %     Sample = parameter();
        %     mu0=Sample.mu;
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

        N_cen = Sample.N_cen;%中心区长度

        H = kron(speye(N_cen), T_0) + kron(diag(ones(N_cen-1,1), 1), T_x) + kron(diag(ones(N_cen-1,1), -1), T_x');
        if periodcity == 1
            H = H + kron( diag(1, N_cen-1), T_x ) + kron( diag(1, -N_cen+1), T_x' );
        end
        GR = inv(speye(4*N_cen)*(EF+1i*eta) - H);
        Gless = -(EF<=0)*(GR - GR');
        % result=diag(kron(eye(N_cen),[0,0,0,1])*Gless);
        result=diag( kron(eye(N_cen),[0,0,0,1])*Gless*kron(eye(N_cen),[1;0;0;0]) );
    end
end