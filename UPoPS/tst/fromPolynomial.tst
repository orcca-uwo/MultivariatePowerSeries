#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));
kernelopts(opaquemodules=false):
UPoPSObject := MultivariatePowerSeries:-UnivariatePolynomialOverPowerSeriesObject:
kernelopts(opaquemodules=true):

# Test (From/To)Polynomial
Try("test 1", UPoPSObject:-Degree(UPoPSObject:-FromPolynomial(0, 'x')), 0);
Try("test 2", UPoPSObject:-ApproximatelyEqual(UPoPSObject:-Zero(), UPoPSObject:-FromPolynomial(0, 'x'), ':-force'), true);

Try("test 3", UPoPSObject:-Degree(UPoPSObject:-FromPolynomial(1, 'x')), 0);
Try("test 4", UPoPSObject:-ApproximatelyEqual(UPoPSObject:-One(), UPoPSObject:-FromPolynomial(1, 'x'), ':-force'), true);

Try[testnoerror]("test 5", UPoPSObject:-FromPolynomial(13*x^10*z + 2*x*y + 2*x + 1, 'x'), 'assign'='upP');
Try[testnoerror]("test 6", UPoPSObject:-Truncate(upP), 'assign'='p');
Try("test 7", evalb(expand(p) = 13*x^10*z + 2*x*y + 2*x + 1), true);

Try[testnoerror]("test 8", UPoPSObject:-Truncate(UPoPSObject:-FromPolynomial(0, 'x')), 'assign'='pp');
Try("test 9", evalb(expand(pp) = 0), true);

Try[testnoerror]("test 10", UPoPSObject:-Truncate(UPoPSObject:-FromPolynomial(-10/3, 'x')), 'assign'='pp');
Try("test 11", evalb(expand(pp) = -10/3), true);

Try("test 12", UPoPSObject:-Degree(upP), 10);
Try("test 13", IsUnit(UPoPSObject:-GetCoefficient(upP, 1)), true);

Try[testnoerror]("test 14", UPoPSObject:-FromPolynomial(2*x*y + 13*z + 1, 'x'), 'assign'='upQ');
Try[testnoerror]("test 15", UPoPSObject:-Truncate(upQ), 'assign'='q');
Try("test 16", evalb(expand(q) = 2*x*y + 13*z + 1), true);

####################################################################

Try("test 17", Degree(UnivariatePolynomialOverPowerSeries(0, 'x')), 0);
Try("test 18", ApproximatelyEqual(UnivariatePolynomialOverPowerSeries(0), UnivariatePolynomialOverPowerSeries(0, 'x'), ':-force'), true);

Try("test 19", Degree(UnivariatePolynomialOverPowerSeries(1, 'x')), 0);
Try("test 20", ApproximatelyEqual(UnivariatePolynomialOverPowerSeries(1), UnivariatePolynomialOverPowerSeries(1, 'x'), ':-force'), true);

Try[testnoerror]("test 21", UnivariatePolynomialOverPowerSeries(13*x^10*z + 2*x*y + 2*x + 1, 'x'), 'assign'='upP');
Try[testnoerror]("test 22", Truncate(upP), 'assign'='p');
Try("test 23", evalb(expand(p) = 13*x^10*z + 2*x*y + 2*x + 1), true);

Try[testnoerror]("test 24", Truncate(UnivariatePolynomialOverPowerSeries(0, 'x')), 'assign'='pp');
Try("test 25", evalb(expand(pp) = 0), true);

Try[testnoerror]("test 26", Truncate(UnivariatePolynomialOverPowerSeries(-10/3, 'x')), 'assign'='pp');
Try("test 27", evalb(expand(pp) = -10/3), true);

Try("test 28", Degree(upP), 10);
Try("test 29", IsUnit(GetCoefficient(upP, 1)), true);

Try[testnoerror]("test 30", UnivariatePolynomialOverPowerSeries(2*x*y + 13*z + 1, 'x'), 'assign'='upQ');
Try[testnoerror]("test 31", Truncate(upQ), 'assign'='q');
Try("test 32", evalb(expand(q) = 2*x*y + 13*z + 1), true);

#end test
