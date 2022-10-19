#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));
kernelopts(opaquemodules=false):
UPoPSObject := MultivariatePowerSeries:-UnivariatePolynomialOverPowerSeriesObject:
kernelopts(opaquemodules=true):

#Setup
Try[testnoerror]("test 1", UPoPSObject:-FromPolynomial(13*x^10*z + 2*x*y + 2*x + 1, 'x'), 'assign'='upP');
Try[testnoerror]("test 2", UPoPSObject:-Truncate(upP), 'assign'='p');
Try[testnoerror]("test 3", UPoPSObject:-FromPolynomial(2*x*y + 13*z + 1, 'x'), 'assign'='upQ');
Try[testnoerror]("test 4", UPoPSObject:-Truncate(upQ), 'assign'='q');

# Test Product
Try[testnoerror]("test 5", UPoPSObject:-Zero(), 'assign'='upZ');
Try[testnoerror]("test 5", UPoPSObject:-BinaryMultiply(upP, upZ), 'assign'='upProd');
Try[testnoerror]("test 6", UPoPSObject:-Truncate(upProd), 'assign'='prod');
Try("test 7", evalb(expand(prod) = 0), true);

Try[testnoerror]("test 8", UPoPSObject:-BinaryMultiply(upP, upQ), 'assign'='upProd');
Try[testnoerror]("test 9", UPoPSObject:-Truncate(upProd, 20), 'assign'='prod');
Try("test 10", evalb(expand(prod) = expand(p*q)), true);

Try[testnoerror]("test 11", UPoPSObject:-BinaryMultiply(upP, UPoPSObject:-One()), 'assign'='upProd');
Try[testnoerror]("test 12", UPoPSObject:-Truncate(upProd, 10), 'assign'='prod');
Try("test 13", evalb(prod = p), true);

####################################################################

#Setup
Try[testnoerror]("test 14", UnivariatePolynomialOverPowerSeries(13*x^10*z + 2*x*y + 2*x + 1, 'x'), 'assign'='upP');
Try[testnoerror]("test 15", Truncate(upP), 'assign'='p');
Try[testnoerror]("test 16", UnivariatePolynomialOverPowerSeries(2*x*y + 13*z + 1, 'x'), 'assign'='upQ');
Try[testnoerror]("test 17", Truncate(upQ), 'assign'='q');

# Test Product
Try[testnoerror]("test 18", Multiply(upP, UnivariatePolynomialOverPowerSeries(0)), 'assign'='upProd');
Try[testnoerror]("test 19", Truncate(upProd), 'assign'='prod');
Try("test 20", evalb(expand(prod) = 0), true);

Try[testnoerror]("test 21", Multiply(upP, upQ), 'assign'='upProd');
Try[testnoerror]("test 22", Truncate(upProd, 20), 'assign'='prod');
Try("test 23", evalb(expand(prod) = expand(p*q)), true);

Try[testnoerror]("test 24", Multiply(upP, UnivariatePolynomialOverPowerSeries(1)), 'assign'='upProd');
Try[testnoerror]("test 25", Truncate(upProd, 10), 'assign'='prod');
Try("test 26", evalb(prod = p), true);

Try[testnoerror]("test 27", MultivariatePowerSeries:-UnivariatePolynomialOverPowerSeries, 'assign'='UPoPS');
Try("test 28", ApproximatelyEqual(Multiply(UPoPS(1), UPoPS(2), UPoPS(3), UPoPS(4), UPoPS(x, x), UPoPS(y, x)), UPoPS(24*x*y, x), 20, ':-force'), true);
Try("test 29", ApproximatelyEqual(Multiply(UPoPS(1), UPoPS(1)), UPoPS(1), 1, ':-force'), true);
Try("test 30", ApproximatelyEqual(Multiply(UPoPS(1), UPoPS(2)), UPoPS(2), 1, ':-force'), true);

#end test
