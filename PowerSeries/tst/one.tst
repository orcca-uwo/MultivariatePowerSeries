#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PSO := MultivariatePowerSeries:-PowerSeriesObject:
kernelopts(opaquemodules=true):

# Test One
Try[testnoerror]("test 1", PSO:-One(), 'assign'='psOne');
Try("test 2", PSO:-ApproximatelyZero(psOne), false);
Try("test 3", PSO:-IsUnit(psOne), true);
Try("test 4", PSO:-ApproximatelyEqual(psOne, PSO:-Zero()), false);
Try("test 5", PSO:-ApproximatelyEqual(psOne, psOne), true);
Try("test 6", PSO:-HomogeneousPart(psOne, 0), 1);
Try("test 7", PSO:-HomogeneousPart(psOne, 1), 0);
Try("test 8", PSO:-HomogeneousPart(psOne, 10), 0);

Try("test 9",  PSO:-GetCoefficient(psOne, 1), 1);
Try("test 10", PSO:-GetCoefficient(psOne, x), 0);
Try("test 11", PSO:-ApproximatelyEqual(psOne, PSO:-Zero(), 0), false);
Try("test 12", PSO:-ApproximatelyEqual(psOne, PSO:-Zero(), -1), true);
Try("test 13", PSO:-ApproximatelyEqual(psOne, PSO:-Zero(), -2), false);


####################################################################

Try[testnoerror]("test 14", PowerSeries(1), 'assign'='psOne');
Try("test 15", ApproximatelyZero(psOne), false);
Try("test 16", IsUnit(psOne), true);
Try("test 17", ApproximatelyEqual(psOne, PowerSeries(0)), false);
Try("test 18", ApproximatelyEqual(psOne, psOne), true);
Try("test 19", HomogeneousPart(psOne, 0), 1);
Try("test 20", HomogeneousPart(psOne, 1), 0);
Try("test 21", HomogeneousPart(psOne, 10), 0);

Try("test 22", GetCoefficient(psOne, 1), 1);
Try("test 23", GetCoefficient(psOne, x), 0);
Try("test 24", ApproximatelyEqual(psOne, PowerSeries(0), 0), false);
Try("test 25", ApproximatelyEqual(psOne, PowerSeries(0), -1), true);
Try("test 26", ApproximatelyEqual(psOne, PowerSeries(0), -2), false);
Try("test 27", ApproximatelyEqual( Multiply(PowerSeries(1), SumOfAllMonomials([x])), SumOfAllMonomials([x]), 10), true);
Try("test 28", ApproximatelyEqual( Add(PowerSeries(1), PowerSeries(1)), PowerSeries(2)), true);

#end test