#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PSO := MultivariatePowerSeries:-PowerSeriesObject:
kernelopts(opaquemodules=true):

# Test Identity
Try[testnoerror]("test 1", PSO:-Identity(x), 'assign'='psI');
Try("test 2", PSO:-ApproximatelyZero(psI), false);
Try("test 3", PSO:-IsUnit(psI), false);
Try("test 4", PSO:-ApproximatelyEqual(psI, PSO:-Zero(), 0), true);
Try("test 5", PSO:-ApproximatelyEqual(psI, PSO:-Zero()), false);
Try("test 6", PSO:-ApproximatelyEqual(psI, PSO:-One()), false);
Try("test 7", PSO:-ApproximatelyEqual(psI, psI), true);
Try("test 8", PSO:-HomogeneousPart(psI, 0), 0);
Try("test 9", PSO:-HomogeneousPart(psI, 1), x);
Try("test 10", PSO:-HomogeneousPart(psI, 2), 0);
Try("test 11", PSO:-GetCoefficient(psI, 1), 0);
Try("test 12", PSO:-GetCoefficient(psI, x), 1);
Try("test 13", PSO:-GetCoefficient(psI, x*y*z), 0);


####################################################################

Try[testnoerror]("test 14", PowerSeries(x), 'assign'='psI');
Try("test 15", ApproximatelyZero(psI), false);
Try("test 16", IsUnit(psI), false);
Try("test 17", ApproximatelyEqual(psI, PowerSeries(0), 0), true);
Try("test 18", ApproximatelyEqual(psI, PowerSeries(0)), false);
Try("test 19", ApproximatelyEqual(psI, PowerSeries(1)), false);
Try("test 20", ApproximatelyEqual(psI, psI), true);
Try("test 21", HomogeneousPart(psI, 0), 0);
Try("test 22", HomogeneousPart(psI, 1), x);
Try("test 23", HomogeneousPart(psI, 2), 0);
Try("test 24", GetCoefficient(psI, 1), 0);
Try("test 25", GetCoefficient(psI, x), 1);
Try("test 26", GetCoefficient(psI, x*y*z), 0);

#end test