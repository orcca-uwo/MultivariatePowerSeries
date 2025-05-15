#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PSO := MultivariatePowerSeries:-PowerSeriesObject:
from_linear_coefficient_ode := MultivariatePowerSeries:-PowerSeriesObject:-from_linear_coefficient_ode:
kernelopts(opaquemodules=true):
ind_var := t;

# Test ode 1
Try[testnoerror]("test 1", from_linear_coefficient_ode(t, x,{diff(x(t),t)- x(t)= 0,x(0)= 1}), 'assign'='ps');
Try("test 2", HomogeneousPart(ps, 0),1);
Try("test 3", HomogeneousPart(ps, 1), 1*ind_var);
Try("test 4", GetCoefficient(ps,(ind_var)^0), 1);
Try("test 5", GetCoefficient(ps,(ind_var)^3), 1/6);
Try("test 6", HomogeneousPart(ps,2),(ind_var^2)/2);
Try("test 7", HomogeneousPart(ps,9),(ind_var^9)/9!);

# Test ode 2
Try[testnoerror]("test 8", from_linear_coefficient_ode(t, x, {diff(x(t),t,t)+x(t)=0,x(0)=0,D(x)(0)=1}), 'assign'='ps2');
Try("test 9", HomogeneousPart(ps2, 0), 0);
Try("test 10", HomogeneousPart(ps2, 1), 1*ind_var);
Try("test 11", GetCoefficient(ps2,(ind_var)^0), 0);
Try("test 12", GetCoefficient(ps2,(ind_var)^3), -1/6);


#Test ode 3
Try[testnoerror]("test 13",from_linear_coefficient_ode(t, x, {diff(x(t),t,t)+ 6*diff(x(t),t)= -5*x(t),x(0)=3,D(x)(0)=-1}), 'assign'='ps3');
Try("test 14", HomogeneousPart(ps3,0),3);
Try("test 15", HomogeneousPart(ps3,2),-9/2*ind_var^2);

#Test ode 4
Try[testnoerror]("test 16",from_linear_coefficient_ode(t, x, {diff(x(t),t,t)+ 6*diff(x(t),t)+9*x(t)=0,x(0)=3,D(x)(0)=-1}), 'assign'='ps4');
Try("test 17", HomogeneousPart(ps4,0),3);
Try("test 18", HomogeneousPart(ps4,2),-21/2*ind_var^2);

#Test ode 5
Try[testnoerror]("test 19",from_linear_coefficient_ode(t, x, {diff(x(t),t,t,t)-diff(x(t),t,t)-20*x(t)=0,x(0)=0,D(x)(0)=1,(D@@2)(x)(0)=1}), 'assign'='ps5');
Try("test 20",HomogeneousPart(ps5,0),0);
Try("test 21",HomogeneousPart(ps5,1),t);
Try("test 22",HomogeneousPart(ps5,2),1/2*ind_var^2);
Try("test 23", HomogeneousPart(ps5,5),41/120*t^5);


#Test ode 6
Try[testnoerror]("test 24",from_linear_coefficient_ode(t, x, {diff(x(t),t,t,t)-3*diff(x(t),t,t)-diff(x(t),t)+3*x(t)=0,x(0)=1,D(x)(0)=2,(D@@2)(x)(0)=3}), 'assign'='ps6');

#Test ode 7
Try[testnoerror]("test 25",from_linear_coefficient_ode(t, x, {diff(x(t),t,t,t,t)-4*diff(x(t),t,t,t)+8*diff(x(t),t,t)-8*diff(x(t),t)+4*x(t)=0,x(0)=0,D(x)(0)=1,(D@@2)(x)(0)=3}),'assign'='ps7');

#Test ode 8
Try[testnoerror]("test 26", from_linear_coefficient_ode(t, x, {diff(x(t),t,t)-2*x(t)-5 = 2*t,x(0)=1,D(x)(0)=2}),'assign'= 'ps8');

#end test
