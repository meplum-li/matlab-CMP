function [mu fval mu_all] = find_mu()
tic
options = optimset('PlotFcns',{@optimplotx,@optimplotfval},'Display','iter','OutputFcn', @myoutput);
% options = optimset('PlotFcns',{@optimplotx,@optimplotfval},'Display','iter');
mu0 = [0.6,0.7];
mu_all = [];
% [mu fval exitflag output] = fzero(@(mu0) transmission_to_current(mu0),mu0,options);
[mu fval] = fzero(@(mu0) transmission_to_current(mu0),mu0,options);
toc

    function stop = myoutput(x,optimvalues,state)
        stop = false;
        if isequal(state,'iter')
            mu_all = [mu_all; x];
        end
    end
end
