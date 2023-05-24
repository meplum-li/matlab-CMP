clearvars
Sample = parameter();
%%%能带结构
mu=Sample.mu;
h = Sample.h;
delta = Sample.delta;
alphaR = Sample.alphaR;
sigma0 = eye(2);
sigmaX=[0,1;1,0];
sigmaY=[0,-1i;1i,0];
sigmaZ=[1,0;0,-1];
T_0 = (2-mu)*kron(sigmaZ, eye(2)) + h*kron(sigmaZ,sigmaZ) - delta* kron(sigmaY, sigmaY);
T_x = -1*kron(sigmaZ, eye(2)) + alphaR/(2i)*kron(sigmaZ, sigmaY);
Kx = linspace(-pi,pi,801);
a=1;
ENG = zeros(4, length(Kx));
Sx = zeros(4,length(Kx));
Sy = zeros(4,length(Kx));
Sz = zeros(4,length(Kx));
for kn = 1 : length(Kx)
    HTT = T_0 + T_x * exp( 1i * Kx(kn) * a) +  T_x' * exp( - 1i * Kx(kn) * a);
   [V,D] =eig(HTT/2+HTT'/2);
   Sx(:,kn) = diag(V'*kron(sigma0,sigmaX)*V);
   Sy(:,kn) = diag(V'*kron(sigma0,sigmaY)*V);
   Sz(:,kn) = diag(V'*kron(sigma0,sigmaZ)*V);
   ENG(:,kn) =diag(D);
end

f=figure;
f.Position(3:4) = [1800 600];
colormap winter
% plot(Kx, ENG,'-','Color', 'k','LineWidth',1.5)
ENG(:,end)=NaN;
subplot(1,3,1)
for ii = 1 : 4
patch(Kx,ENG(ii,:),Sx(ii,:),'EdgeColor','interp','LineWidth',3)
hold on
end
set(gca, 'FontSize', 20);
colorbar
clim([-1,1])
subplot(1,3,2)
for ii = 1 : 4
patch(Kx,ENG(ii,:),Sy(ii,:),'EdgeColor','interp','LineWidth',3)
hold on
end
set(gca, 'FontSize', 20);
colorbar
clim([-1,1])
subplot(1,3,3)
for ii = 1 : 4
patch(Kx,ENG(ii,:),Sz(ii,:),'EdgeColor','interp','LineWidth',3)
hold on
end
set(gca, 'FontSize', 20);
colorbar
clim([-1,1])
%%
%%%能谱结构
N_cen = 100;
HH = kron(eye(N_cen), T_0) + kron(diag(ones(N_cen-1,1), 1), T_x) + kron(diag(ones(N_cen-1,1), -1), T_x');
ENG = eig(HH);
figure(Position=[100,1000,1800,500])
subplot(1,3,1)
plot(ENG,'s','Color', 'k','LineWidth',1.5)
set(gca, 'FontSize', 20);
ylim([-0.3,0.3])

subplot(1,3,2)
ylim([-0.3,0.3])
[V,~]=eigs(HH,2,'smallestabs');
plot(kron(eye(N_cen),[1,1,1,1])*(conj(V(:,1)).*V(:,1)),'Color', 'k','LineWidth',3)
set(gca, 'FontSize', 20);

subplot(1,3,3)
ylim([-0.3,0.3])
[V,~]=eigs(HH,2,'smallestabs');
plot(kron(eye(N_cen),[1,1,1,1])*(conj(V(:,2)).*V(:,2)),'Color', 'k','LineWidth',3)
set(gca, 'FontSize', 20);
%%%边缘态的分布
