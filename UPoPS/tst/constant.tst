#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));
kernelopts(opaquemodules=false):
UPoPSObject := MultivariatePowerSeries:-UnivariatePolynomialOverPowerSeriesObject:
kernelopts(opaquemodules=true):

# Test Constant 
Try[testnoerror]("test 1", UPoPSObject:-Constant(5), 'assign'='upConst');
Try("test 2", UPoPSObject:-Degree(upConst), 0);
Try("test 3", IsUnit(UPoPSObject:-GetCoefficient(upConst, 0)), true);
Try[testerror]("test 4", UPoPSObject:-GetCoefficient(upConst, 1), "invalid input: degree out of range");
Try("test 5", UPoPSObject:-ApproximatelyEqual(upConst, UPoPSObject:-Zero(), ':-force'), false);
Try("test 6", UPoPSObject:-ApproximatelyEqual(upConst, UPoPSObject:-One(), ':-force'), false);
Try("test 7", UPoPSObject:-ApproximatelyEqual(upConst, upConst, ':-force'), true);

####################################################################

Try[testnoerror]("test 8",  UnivariatePolynomialOverPowerSeries(5), 'assign'='upConst');
Try("test 9", Degree(upConst), 0);
Try("test 10", IsUnit(GetCoefficient(upConst, 0)), true);
Try[testerror]("test 11", GetCoefficient(upConst, 1), "invalid input: degree out of range");
Try("test 12", ApproximatelyEqual(upConst, UnivariatePolynomialOverPowerSeries(0), ':-force'), false);
Try("test 13", ApproximatelyEqual(upConst, UnivariatePolynomialOverPowerSeries(1), ':-force'), false);
Try("test 14", ApproximatelyEqual(upConst, upConst, ':-force'), true);

#end test
