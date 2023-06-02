clearvars
iter_max = 150;
RelTol = 1e-5;
AbsTol = 1e-10;
delta_all = zeros(1,iter_max+1);
RelErr = zeros(2,iter_max+1);
AbsErr = zeros(2,iter_max+1);
mu0 = 0.6;
Ui = 1.74912443;
% Ui = 1.7480;
delta0 = 0.1529;
delta_curr = delta0;
ElapsedTime = toc;
info = ['iteration','   current delta','    AbsErr','   RelErr','   Time Cost\n'];
fprintf('%%------------------------\n')
fprintf(info)
fprintf('%8d %13.5f %16.E %7.E %9.2f\n', [0,delta0,NaN, NaN,NaN])
for ii = 1 : iter_max
    tic
    % fixed point estimate at first iteration
    delta_all(ii) = delta_curr;
    delta_next = zero_bias_delta(Ui, mu0, abs(delta_curr));
    RelErr(:,ii) = [abs(real(delta_next - delta_curr))/real(delta_next);abs(imag(delta_next - delta_curr))/imag(delta_next)];
    AbsErr(:,ii) = [abs(real(delta_next - delta_curr));abs(imag(delta_next - delta_curr))];
    if ( (RelErr(1,ii)<RelTol)||(AbsErr(1,ii)<AbsTol) )&&( (RelErr(2,ii)<RelTol)||(AbsErr(2,ii)<AbsTol) )
        ElapsedTime = toc;
        info = ['iteration','   current delta','    AbsErr','   RelErr','   Time Cost\n'];
        fprintf('%%------------------------\n')
        fprintf(info)
        fprintf('%8d %13.5f %16.E %7.E %9.2f\n', [ii,delta_next,AbsErr(2,ii),RelErr(1,ii),ElapsedTime])
        break
    end
    ElapsedTime = toc;
    info = ['iteration','   current delta','    AbsErr','   RelErr','   Time Cost\n'];
    fprintf('%%------------------------\n')
    fprintf(info)
    fprintf('%8d %13.5f %16.E %7.E %9.2f\n', [ii,delta_next,AbsErr(2,ii),RelErr(1,ii),ElapsedTime])
    delta_curr = delta_next;
end
% converged fixed point
delta = delta_next;
delta_all(ii+1) = delta;
delta_all = delta_all(1:(ii+1));
RelErr = RelErr(:,1:ii);
AbsErr = AbsErr(:,1:ii);
