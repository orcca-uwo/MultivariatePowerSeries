#test
with(TestTools):
Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PSO := MultivariatePowerSeries:-PowerSeriesObject:
kernelopts(opaquemodules=true):

# Test Zero
Try[testnoerror]("test 1", PSO:-Zero(), 'assign'='psZero');
Try("test 2", PSO:-ApproximatelyZero(psZero), true);
Try("test 3", PSO:-IsUnit(psZero), false);
Try("test 4", PSO:-ApproximatelyEqual(psZero, psZero), true);
Try("test 5", PSO:-HomogeneousPart(psZero, 0), 0);
Try("test 6", PSO:-HomogeneousPart(psZero, 1), 0);
Try("test 7", PSO:-HomogeneousPart(psZero, 10), 0);
Try("test 8", PSO:-GetCoefficient(psZero, y^1), 0);

####################################################################

Try[testnoerror]("test 9", PowerSeries(0), 'assign'='psZero');
Try("test 10", ApproximatelyZero(psZero), true);
Try("test 11", IsUnit(psZero), false);
Try("test 12", ApproximatelyEqual(psZero, psZero), true);
Try("test 13", HomogeneousPart(psZero, 0), 0);
Try("test 14", HomogeneousPart(psZero, 1), 0);
Try("test 15", HomogeneousPart(psZero, 10), 0);
Try("test 16", GetCoefficient(psZero, y^1), 0);
Try("test 17", ApproximatelyEqual( Subtract(PowerSeries(1 + x^2), PowerSeries(1 + x^2)), psZero), true);
Try("test 18", ApproximatelyEqual( Add(PowerSeries(0), PowerSeries(0), PowerSeries(0)), psZero), true);
Try("test 19", ApproximatelyEqual( Multiply(PowerSeries(0), SumOfAllMonomials([x]), GeometricSeries([x, y])), psZero), true);

#end test 