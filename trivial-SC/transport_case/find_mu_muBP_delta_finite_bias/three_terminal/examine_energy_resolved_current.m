clearvars
Sample = parameter();
current_RelTol = Sample.current.RelTol;
current_AbsTol = Sample.current.AbsTol;
EF = -0.04;
mu0 = 0.65;
[i_vec] = i_E(Sample, EF, mu0);
% ub = max(abs(Sample.A_mu_exU-mu0),abs(Sample.A_mu_exD-mu0)); 
% i_vec=integral(@(EF) i_E(Sample, EF, mu0),-ub,ub,"ArrayValued",true,'RelTol',current_RelTol,'AbsTol',current_AbsTol);
fprintf('%6.3E  %6.3E  %6.3E\n',i_vec(1:3))

%%
function [i_vec] = i_E(Sample, EF, mu0)
        h = Sample.h;
        delta = Sample.delta;
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
        A_mu_exBP = Sample.A_mu_exBP;
        mu_exU = A_mu_exU - mu0;%relative value
        mu_exD = A_mu_exD - mu0;
        mu_exBP = A_mu_exBP - mu0;
        gammaU = Sample.gammaU;%上导线的gamma
        gammaD = Sample.gammaD;%上导线的gamma
        gammaBP = Sample.gammaBP;

%         TT = zeros(6,6);%不同能量对应的透射系数


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
        grlead_lead= inv(eye(4*N_cen)*EF - H-SigmaU - SigmaD - SigmaBP);
        %%% 计算transmission 矩阵
        GAMMA_T = [GammaU; GammaD; GammaBP];%for calculating transmission
        temT = GAMMA_T*grlead_lead*[GammaU, GammaD, GammaBP].*conj( kron(ones(3),grlead_lead) );%稀疏矩阵的写法提速10倍，非稀疏的就不用写成稀疏矩阵了
        % %保留自旋自由度
        % sumM=blkdiag(eye(8));
        % T = real(sumM*temT*sumM');
        %自旋缩并
        sumM=kron(speye(2*3*N_cen), ones(1,2));
        temTT = real(sumM*temT*sumM');
        %缩并电极内部指标，只剩下up、down、BP以及各自的电子空穴指标
        sumM=kron(  speye(3), kron( ones(1, N_cen), speye(2) )  );
        TT= real(sumM*temTT*sumM');
        fUe = (EF<mu_exU);
        fUh = (EF<-mu_exU);
        fDe = (EF<mu_exD);
        fDh = (EF<-mu_exD);
        fBPe = (EF<mu_exBP);
        fBPh = (EF<-mu_exBP);
        i_vec(1) = sum( (fUe - [fUe,fUh, fDe, fDh, fBPe, fBPh]) *TT(:,1) ) ...
            -sum( (fUh - [fUe, fUh, fDe, fDh, fBPe, fBPh]) *TT(:,2) );%U
        i_vec(2) = sum( (fDe - [fUe, fUh, fDe, fDh, fBPe, fBPh]) *TT(:,3) ) ...
            -sum( (fDh - [fUe, fUh, fDe, fDh, fBPe, fBPh]) *TT(:,4) );%D
        i_vec(3) = sum( (fBPe - [fUe, fUh, fDe, fDh, fBPe, fBPh]) *TT(:,5) ) ...
            -sum( (fBPh - [fUe, fUh, fDe, fDh, fBPe, fBPh]) *TT(:,6) );%BP
    end