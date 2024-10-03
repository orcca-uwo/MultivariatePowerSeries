#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
ps1O := MultivariatePowerSeries:-PowerSeriesObject:
from_polynom_coefficient_ode := MultivariatePowerSeries:-PowerSeriesObject:-from_polynom_coefficient_ode:
ode := MultivariatePowerSeries:-OdeManipulator:-OdeObject;
kernelopts(opaquemodules=true):


#Test on ode1
Try[testnoerror]("test 1", from_polynom_coefficient_ode(t, x,{diff(x(t),t)- 2*t*x(t)= 0, x(0)=1}), 'assign'='ps1');
Try("test 2",HomogeneousPart(ps1,0),1);
Try("test 3",HomogeneousPart(ps1,1),0);
Try("test 4",HomogeneousPart(ps1,2),t^2);
Try("test 5",HomogeneousPart(ps1,3),0);
Try("test 6",HomogeneousPart(ps1,4),1/2*t^4);

#Test on ode2
Try[testnoerror]("test 7", from_polynom_coefficient_ode(t,x,{diff(x(t),t,t)-t^2*x(t)=0}),'assign'='ps2');
Try("test 8", HomogeneousPart(ps2,0),x(0));
Try("test 9",HomogeneousPart(ps2,5),1/20*D(x)(0)*t^5);


#Test on ode3
Try[testnoerror]("test 10", from_polynom_coefficient_ode(t,x,{diff(x(t),t,t)+ 3*t*diff(x(t),t)- x(t)= 0,x(0)= 2, D(x)(0)= 0}),'assign'='ps3');
Try("test 11", HomogeneousPart(ps3,0),2);
Try("test 12", HomogeneousPart(ps3,3),0);
Try("test 13", HomogeneousPart(ps3,6),11/72*t^6);

#Test on ode4
Try[testerror]("test 14", from_polynom_coefficient_ode(t,x,{t^2*diff(x(t),t,t)+ x(t)= 0}),'assign'='ps4');


Try[testnoerror]("test 15", from_polynom_coefficient_ode(t, x,{diff(x(t),t)- 2*t*x(t)= 0, x(0)=k}), 'assign'='ps');
Try("test 16", HomogeneousPart(ps,0),k);

k := 3:
Try("test 17", HomogeneousPart(ps,0),3);
Try("test 18", HomogeneousPart(ps,2),3*t^2);

#end test