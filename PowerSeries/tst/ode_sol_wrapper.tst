#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
ps1O := MultivariatePowerSeries:-PowerSeriesObject:
ode_sol_wrapper := MultivariatePowerSeries:-PowerSeriesObject:-ode_sol_wrapper:
kernelopts(opaquemodules=true):

# Test ode 1 constant coefficient
Try[testnoerror]("test 1", ode_sol_wrapper(t, x, {diff(x(t),t,t)+x(t)=0,x(0)=0,D(x)(0)=1}), 'assign'='ps2');
Try("test 2", HomogeneousPart(ps2, 0), 0);
Try("test 3", HomogeneousPart(ps2, 1), 1*t);
Try("test 4", GetCoefficient(ps2,(t)^0), 0);
Try("test 5", GetCoefficient(ps2,(t)^3), -1/6);

#Test on ode 2 polynomial coefficient
Try[testnoerror]("test 6", ode_sol_wrapper(t,x,{diff(x(t),t,t)-t^2*x(t)=0}),'assign'='ps2');
Try("test 7", HomogeneousPart(ps2,0),x(0));
Try("test 8",HomogeneousPart(ps2,5),1/20*D(x)(0)*t^5);

#end test