#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
ps1O := MultivariatePowerSeries:-PowerSeriesObject:
solution_around_exp_pnt := MultivariatePowerSeries:-PowerSeriesObject:-solution_around_exp_pnt:
ode := MultivariatePowerSeries:-OdeManipulator:-OdeObject;
kernelopts(opaquemodules=true):

Try[testnoerror]("test 1", solution_around_exp_pnt(t, x,{diff(x(t),t,t)- sin(t)*x(t)= 0, x(1)=exp(1),D(x)(1)=1/2},exp_pnt = 1), 'assign'='ps1');

Try[testnoerror]("test 2",solution_around_exp_pnt(t,x,{diff(x(t),t,t) - sin(t)*x(t) = 0, x(Pi) = 1, D(x)(Pi) = 0},exp_pnt = Pi),'assign'='ps2');


#end test


