#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));
kernelopts(opaquemodules=false):
UPoPSObject := MultivariatePowerSeries:-UnivariatePolynomialOverPowerSeriesObject:
kernelopts(opaquemodules=true):

# Test One
Try[testnoerror]("test 1", UPoPSObject:-One(), 'assign'='upOne');
Try("test 2", UPoPSObject:-Degree(upOne), 0);
Try("test 3", IsUnit(UPoPSObject:-GetCoefficient(upOne, 0)), true);
Try[testerror]("test 4", UPoPSObject:-GetCoefficient(upOne, 1), "invalid input: degree out of range");
Try("test 5", UPoPSObject:-ApproximatelyEqual(upOne, UPoPSObject:-Zero(), ':-force'), false);
Try("test 6", UPoPSObject:-ApproximatelyEqual(upOne, upOne, ':-force'), true);

####################################################################

Try[testnoerror]("test 7", UnivariatePolynomialOverPowerSeries(1), 'assign'='upOne');
Try("test 8", Degree(upOne), 0);
Try("test 9", IsUnit(GetCoefficient(upOne, 0)), true);
Try[testerror]("test 10", GetCoefficient(upOne, 1), "invalid input: degree out of range");
Try("test 11", ApproximatelyEqual(upOne, UnivariatePolynomialOverPowerSeries(0), ':-force'), false);
Try("test 12", ApproximatelyEqual(upOne, upOne, ':-force'), true);

#end test
