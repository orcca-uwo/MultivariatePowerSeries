#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PSO := MultivariatePowerSeries:-PowerSeriesObject:
kernelopts(opaquemodules=true):

# Test Constant
Try[testnoerror]("test 1", PSO:-Constant(0), 'assign'='psConst');
Try("test 2", PSO:-ApproximatelyZero(psConst), true);
Try("test 3", PSO:-ApproximatelyEqual(psConst, PSO:-Zero()), true);
Try("test 4", PSO:-ApproximatelyEqual(psConst, PSO:-One()), false);
Try("test 5", PSO:-ApproximatelyEqual(psConst, PSO:-Identity(x)), false);

Try[testnoerror]("test 6", PSO:-Constant(1), 'assign'='psConst');
Try("test 7", PSO:-ApproximatelyZero(psConst), false);
Try("test 8", PSO:-IsUnit(psConst), true);
Try("test 9", PSO:-ApproximatelyEqual(psConst, PSO:-Zero()), false);
Try("test 10", PSO:-ApproximatelyEqual(psConst, PSO:-One()), true);
Try("test 11", PSO:-ApproximatelyEqual(psConst, PSO:-Identity(x)), false);
Try("test 12", PSO:-ApproximatelyEqual(psConst, psConst), true);

Try[testnoerror]("test 13", PSO:-Constant(5), 'assign'='psConst');
Try("test 14", PSO:-HomogeneousPart(psConst, 0), 5);
Try("test 15", PSO:-HomogeneousPart(psConst, 1), 0);
Try("test 16", PSO:-HomogeneousPart(psConst, 10), 0);
Try("test 17", PSO:-GetCoefficient(psConst, 1), 5);
Try("test 18", PSO:-GetCoefficient(psConst, x), 0);
Try("test 19", PSO:-GetCoefficient(psConst, x^5*y^5*z^5), 0);

Try[testnoerror]("test 20", PSO:-Constant(RootOf(z^2 + 1)), 'assign'='psConst');
Try("test 21", PSO:-HomogeneousPart(psConst, 0), RootOf(z^2 + 1));
Try("test 22", PSO:-HomogeneousPart(psConst, 1), 0);

Try[testnoerror]("test 23", PSO:-Constant(sqrt(3)), 'assign'='pssqrt');
Try("test 24", PSO:-HomogeneousPart(pssqrt, 0), RootOf(_Z^2-3,index = 1));
Try("test 25", PSO:-HomogeneousPart(pssqrt, 1), 0);



####################################################################

Try[testnoerror]("test 26",   PowerSeries(0), 'assign'='psConst');
Try("test 27",  ApproximatelyZero(psConst), true);
Try("test 28",  ApproximatelyEqual(psConst,  PowerSeries(0)), true);
Try("test 29",  ApproximatelyEqual(psConst,  PowerSeries(1)), false);
Try("test 30",  ApproximatelyEqual(psConst,  PowerSeries(x)), false);

Try[testnoerror]("test 31",   PowerSeries(1), 'assign'='psConst');
Try("test 32",  ApproximatelyZero(psConst), false);
Try("test 33",  IsUnit(psConst), true);
Try("test 33",  ApproximatelyEqual(psConst,  PowerSeries(0)), false);
Try("test 34",  ApproximatelyEqual(psConst,  PowerSeries(1)), true);
Try("test 35",  ApproximatelyEqual(psConst,  PowerSeries(x)), false);
Try("test 36",  ApproximatelyEqual(psConst, psConst), true);

Try[testnoerror]("test 37",   PowerSeries(5), 'assign'='psConst');
Try("test 38",  HomogeneousPart(psConst, 0), 5);
Try("test 39",  HomogeneousPart(psConst, 1), 0);
Try("test 40",  HomogeneousPart(psConst, 10), 0);
Try("test 41",  GetCoefficient(psConst, 1), 5);
Try("test 42",  GetCoefficient(psConst, x), 0);
Try("test 43",  GetCoefficient(psConst, x^5*y^5*z^5), 0);

Try[testnoerror]("test 44",   PowerSeries(RootOf(z^2 + 1)), 'assign'='psConst');
Try("test 45",  HomogeneousPart(psConst, 0), RootOf(z^2 + 1));
Try("test 46",  HomogeneousPart(psConst, 1), 0);

Try[testnoerror]("test 47",   PowerSeries(sqrt(3)), 'assign'='pssqrt');
Try("test 48",  HomogeneousPart(pssqrt, 0), RootOf(_Z^2-3,index = 1));
Try("test 49",  HomogeneousPart(pssqrt, 1), 0);

#end test
