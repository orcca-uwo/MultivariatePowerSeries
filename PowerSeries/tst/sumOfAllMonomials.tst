#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PSO := MultivariatePowerSeries:-PowerSeriesObject:
kernelopts(opaquemodules=true):


# Test SumOfAllMonomials
Try[testnoerror]("test 1", PSO:-SumOfAllMonomials([x, y]), 'assign'='psSAM');
Try("test 2", PSO:-ApproximatelyZero(psSAM, 'force'), false);
Try("test 3", PSO:-IsUnit(psSAM), true);
Try("test 4", PSO:-ApproximatelyEqual(psSAM, psSAM, ':-force'), true);
Try("test 5", PSO:-ApproximatelyEqual(psSAM, PSO:-GeometricSeries([x, y]), ':-force'), true);
Try("test 6", PSO:-HomogeneousPart(psSAM, 0), 1);
Try("test 7", PSO:-HomogeneousPart(psSAM, 1), x + y);
Try("test 8", PSO:-HomogeneousPart(psSAM, 2), x^2 + x*y + y^2);
Try("test 9", PSO:-GetCoefficient(psSAM, x*y^4), 1);
Try("test 10", degree(PSO:-HomogeneousPart(psSAM, 10)), 10);

####################################################################

Try[testnoerror]("test 11", SumOfAllMonomials([x, y]), 'assign'='psSAM');
Try("test 12", ApproximatelyZero(psSAM, 'force'), false);
Try("test 13", IsUnit(psSAM), true);
Try("test 14", ApproximatelyEqual(psSAM, psSAM, ':-force'), true);
Try("test 15", ApproximatelyEqual(psSAM, GeometricSeries([x, y]), ':-force'), true);
Try("test 16", HomogeneousPart(psSAM, 0), 1);
Try("test 17", HomogeneousPart(psSAM, 1), x + y);
Try("test 18", HomogeneousPart(psSAM, 2), x^2 + x*y + y^2);
Try("test 19", GetCoefficient(psSAM, x*y^4), 1);
Try("test 20", degree(HomogeneousPart(psSAM, 10)), 10);


#end test
