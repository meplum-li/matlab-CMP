function [F, J] = mu_muBP_Delta_Jacobi(Sample, x)
mu0 = x(1);
mu_BP = x(2);%A_mu_exBP
delta = x(3);
Ui = Sample.Ui;
gap_RelTol = Sample.gap.RelTol;%简单用gap的积分精度来全局控制gap和电流的积分
gap_AbsTol = Sample.gap.AbsTol;
int=integral(@(EF) current_diffDelta_jacobi(Sample, EF, mu0, mu_BP, delta),-inf,inf,"ArrayValued",true,'RelTol',gap_RelTol,'AbsTol',gap_AbsTol);
F(1) = sum(int(1:3));
F(2) = int(3);
F(3) = -Ui/(2*pi*1i)*int(4) - delta;
[TT_Ue, Ue_GRGA41] = quantity_given_EF(Sample, Sample.A_mu_exU-mu0, mu0, mu_BP, delta);
[TT_Uh, Uh_GRGA41] = quantity_given_EF(Sample, -Sample.A_mu_exU+mu0, mu0, mu_BP, delta);
[TT_De, De_GRGA41] = quantity_given_EF(Sample, Sample.A_mu_exD-mu0, mu0, mu_BP, delta);
[TT_Dh, Dh_GRGA41] = quantity_given_EF(Sample, -Sample.A_mu_exD+mu0, mu0, mu_BP, delta);
[TT_BPe, BPe_GRGA41] = quantity_given_EF(Sample, mu_BP-mu0, mu0, mu_BP, delta);
[TT_BPh, BPh_GRGA41] = quantity_given_EF(Sample, -mu_BP+mu0, mu0, mu_BP, delta);

J(1,1) = sum(int(5:7)) +sum(sum( -[TT_Ue(1,:);TT_Uh(2,:);TT_De(3,:);TT_Dh(4,:);TT_BPe(5,:);TT_BPh(6,:)]+kron(ones(3),[1,-1;-1,1]).*[TT_Ue(:,1),TT_Uh(:,2),TT_De(:,3),TT_Dh(:,4),TT_BPe(:,5),TT_BPh(:,6)] ));%DF(1)/Dmu0
J(1,2) = sum([TT_BPe(5,:),TT_BPh(6,:)])-sum(sum( kron(ones(3,1),[1,-1;-1,1]).*[TT_BPe(:,5),TT_BPh(:,6)] ));%DF(1)/Dmu_BP
J(1,3) = sum(int(8:10));%DF(1)/Ddelta
J(2,1) = int(10) +sum( -[TT_BPe(5,:),TT_BPh(6,:)] ) + sum(sum( kron(ones(1,3),[1,-1;-1,1]).*[TT_Ue([5,6],1),TT_Uh([5,6],2),TT_De([5,6],3),TT_Dh([5,6],4),TT_BPe([5,6],5),TT_BPh([5,6],6)] ));%DF(2)/Dmu0
J(2,2) = sum([TT_BPe(5,:),TT_BPh(6,:)])-sum(sum( [1,-1;-1,1].*[TT_BPe([5,6],5),TT_BPh([5,6],6)] ));%DF(2)/Dmu_BP
J(2,3) = int(10);%DF(2)/Ddelta
J(3,1) = -Ui/(2*pi*1i)*int(11) -Ui/(2*pi*1i)*(-1i)*sum( kron(ones(3,1),[1;-1]).*[Ue_GRGA41;Uh_GRGA41;De_GRGA41;Dh_GRGA41;BPe_GRGA41;BPh_GRGA41] );%DF(3)/Dmu0
J(3,2) = -Ui/(2*pi*1i)*1i * sum([1;-1].*[BPe_GRGA41;BPh_GRGA41]);%DF(3)/Dmu_BP
J(3,3) = -Ui/(2*pi*1i)*int(12) - 1;%DF(3)/Ddelta

end