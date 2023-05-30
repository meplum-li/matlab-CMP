
% tic
% delta0 = [0.23,0.35];
% [delta fval exitflag output] = fsolve(@(delta) finite_bias_self_consist_delta(delta,mu),delta0,options)
% toc
clearvars
iter_max = 20;
RelTol = 1e-4;
AbsTol = 1e-10;
delta0 = 0.3;
delta_curr = delta0;
delta_all = zeros(1,iter_max+1);
RelErr = zeros(2,iter_max+1);
AbsErr = zeros(2,iter_max+1);
for ii = 1 : iter_max
    % fixed point estimate at first iteration
    delta_all(ii) = delta_curr;
    delta_next = finite_bias_self_consist_delta(abs(delta_curr), 0.6);
    RelErr(:,ii) = [abs(real(delta_next - delta_curr))/real(delta_next);abs(imag(delta_next - delta_curr))/imag(delta_next)];
    AbsErr(:,ii) = [abs(real(delta_next - delta_curr));abs(imag(delta_next - delta_curr))];
    if ( (RelErr(1,ii)<RelTol)||(AbsErr(1,ii)<AbsTol) )&&( (RelErr(2,ii)<RelTol)||(AbsErr(2,ii)<AbsTol) )
        break
    end
    delta_curr = delta_next;
end
    % converged fixed point
    delta = delta_next;
            delta_all(ii+1) = delta;
        delta_all = delta_all(1:(ii+1));
        RelErr = RelErr(:,1:ii);
        AbsErr = AbsErr(:,1:ii);
