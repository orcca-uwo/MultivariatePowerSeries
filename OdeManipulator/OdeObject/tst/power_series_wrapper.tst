#test
##############
## ODEtests ##
##############
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
ODEO := MultivariatePowerSeries:-OdeManipulator:-OdeObject:
PowerSeriesWrapper := MultivariatePowerSeries:-OdeManipulator:-PowerSeriesWrapper:
DummyConstructor := MultivariatePowerSeries:-OdeManipulator:-OdeObject:-DummyConstructor;
kernelopts(opaquemodules=true):

Try[testnoerror]("test 1", Object(PowerSeriesWrapper, PowerSeries(t^2)), 'assign'='p1');
Try[testnoerror]("test 2", Object(PowerSeriesWrapper, PowerSeries(sin(t))), 'assign'='p2');
Try[testnoerror]("test 3", Object(PowerSeriesWrapper, PowerSeries(cos(t))), 'assign'='p3');
Try[testnoerror]("test 3", Object(PowerSeriesWrapper, PowerSeries((sin(t))^2 )), 'assign'='p4');

Try[testnoerror]("test 4", diff(p1, x), 'assign'='dp');
Try("test 5", dp:-Truncate(dp, 5), 0);

Try[testnoerror]("test 6", diff(p1, t), 'assign'='dp');
Try("test 7", dp:-Truncate(dp, 5), 2*t);

Try[testnoerror]("test 8", diff(p2, t), 'assign'='dp');
Try("test 9", dp:-Truncate(dp, 10), p3:-Truncate(p3, 10));

Try[testnoerror]("test 10", diff(p3, t), 'assign'='dp');
Try("test 11", dp:-Truncate(dp, 10), -p2:-Truncate(p2, 10));

Try[testnoerror]("test 12", diff(p4, t), 'assign'='dp');
r :=  2*PowerSeries(sin(t))*PowerSeries(cos(t)):
Try("test 13", dp:-Truncate(dp, 10), Truncate(r,10));

#end test
