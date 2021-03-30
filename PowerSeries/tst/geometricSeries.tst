#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PSO := MultivariatePowerSeries:-PowerSeriesObject:
kernelopts(opaquemodules=true):

# Test GeometricSeries
Try[testnoerror]("test 1", PSO:-GeometricSeries([x, y]), 'assign'='psGeo');
Try("test 2", PSO:-ApproximatelyZero(psGeo), false);
Try("test 3", PSO:-IsUnit(psGeo), true);
Try("test 4", PSO:-ApproximatelyEqual(psGeo, psGeo), true);
Try("test 5", PSO:-HomogeneousPart(psGeo, 0), 1);
Try("test 6", PSO:-HomogeneousPart(psGeo, 1), x + y);
Try("test 7", PSO:-HomogeneousPart(psGeo, 2), x^2 + 2*x*y + y^2);
Try("test 8", PSO:-GetCoefficient(psGeo, x^3*y^2), 10);
Try("test 9", degree(PSO:-HomogeneousPart(psGeo, 10)), 10);


####################################################################

# Test GeometricSeries
Try[testnoerror]("test 10", GeometricSeries([x, y]), 'assign'='psGeo');
Try("test 11", ApproximatelyZero(psGeo), false);
Try("test 12", IsUnit(psGeo), true);
Try("test 13", ApproximatelyEqual(psGeo, psGeo), true);
Try("test 14", HomogeneousPart(psGeo, 0), 1);
Try("test 15", HomogeneousPart(psGeo, 1), x + y);
Try("test 16", HomogeneousPart(psGeo, 2), x^2 + 2*x*y + y^2);
Try("test 17", GetCoefficient(psGeo, x^3*y^2), 10);
Try("test 18", degree(HomogeneousPart(psGeo, 10)), 10);

#end test