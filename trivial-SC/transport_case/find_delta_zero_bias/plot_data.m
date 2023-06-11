figure
plot(gamma,delta,'*-','LineWidth',2)
fun = @(x) real(x*sqrt(1-2*gamma./x)-delta);
x0 = 0.33;
x = lsqnonlin(fun,x0);
x=0.3;
plot(gamma,delta,'ko',gamma,x*sqrt(1-2*gamma./x),'b-')
legend('Data','Best fit')
xlabel('gamma')
ylabel('exp(-tx)')
