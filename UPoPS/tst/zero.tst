#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));
kernelopts(opaquemodules=false):
UPoPSObject := MultivariatePowerSeries:-UnivariatePolynomialOverPowerSeriesObject:
kernelopts(opaquemodules=true):

# Test Zero
Try[testnoerror]("test 1", UPoPSObject:-Zero(), 'assign'='upZero');
Try("test 2", UPoPSObject:-Degree(upZero), 0);
Try("test 3", UPoPSObject:-Truncate(upZero, 10), 0);
Try[testerror]("test 4", UPoPSObject:-GetCoefficient(upZero, 1), "invalid input: degree out of range");
Try("test 5", UPoPSObject:-ApproximatelyEqual(upZero, upZero), true);

####################################################################

Try[testnoerror]("test 6", UnivariatePolynomialOverPowerSeries(0), 'assign'='upZero');
Try("test 7", Degree(upZero), 0);
Try("test 8", ApproximatelyZero(GetCoefficient(upZero, 0)), true);
Try[testerror]("test 9", GetCoefficient(upZero, 1), "invalid input: degree out of range");
Try("test 10", ApproximatelyEqual(upZero, upZero), true);

#end test