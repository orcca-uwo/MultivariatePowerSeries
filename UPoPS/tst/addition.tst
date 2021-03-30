#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));
kernelopts(opaquemodules=false):
UPoPSObject := MultivariatePowerSeries:-UnivariatePolynomialOverPowerSeriesObject:
kernelopts(opaquemodules=true):

#Setup 
Try[testnoerror]("test 1", UPoPSObject:-One(), 'assign'='upOne');
Try[testnoerror]("test 2", UPoPSObject:-FromPolynomial(13*x^10*z + 2*x*y + 2*x + 1, 'x'), 'assign'='upP');
Try[testnoerror]("test 3", UPoPSObject:-Truncate(upP), 'assign'='p');
Try[testnoerror]("test 4", UPoPSObject:-FromPolynomial(2*x*y + 13*z + 1, 'x'), 'assign'='upQ');
Try[testnoerror]("test 5", UPoPSObject:-Truncate(upQ), 'assign'='q');

# Test Add
Try[testnoerror]("test 6", UPoPSObject:-BinaryAdd(upP, upQ), 'assign'='upPaQ');
Try[testnoerror]("test 7", UPoPSObject:-Truncate(upPaQ, 10), 'assign'='paq');
Try("test 8", evalb(expand(paq) = expand(p+q)), true);

Try[testnoerror]("test 9", UPoPSObject:-BinaryAdd(upP, upP), 'assign'='upPaP');
Try[testnoerror]("test 10", UPoPSObject:-Truncate(upPaP, 10), 'assign'='pap');
Try("test 11", evalb(expand(pap) = expand(2*p)), true);

Try[testnoerror]("test 12", UPoPSObject:-BinaryAdd(upPaP, upPaQ), 'assign'='upPaPaPaQ');
Try[testnoerror]("test 13", UPoPSObject:-Truncate(upPaPaPaQ, 10), 'assign'='papapaq');
Try("test 14", evalb(expand(papapaq) = expand(3*p + q)), true);

# Test Sub
Try[testnoerror]("test 15", UPoPSObject:-BinarySub(upPaPaPaQ, upP), 'assign'='upPaPaQ');
Try[testnoerror]("test 16", UPoPSObject:-BinarySub(upPaPaQ, upP), 'assign'='upPaQ');
Try[testnoerror]("test 17", UPoPSObject:-BinarySub(upPaQ, upP), 'assign'='upaQ');
Try[testnoerror]("test 18", UPoPSObject:-BinarySub(upaQ, upQ), 'assign'='upsub');
Try[testnoerror]("test 19", UPoPSObject:-Truncate(upsub, 10), 'assign'='sub');
Try("test 20", evalb(expand(sub) = 0), true);

# Test Negate 
Try[testnoerror]("test 21", UPoPSObject:-BinaryAdd(upP, UPoPSObject:-Negate(upP)), 'assign'='upPmP');
Try[testnoerror]("test 22", UPoPSObject:-Truncate(upPmP), 'assign'='sub');
Try("test 23", evalb(expand(sub) = 0), true);

Try[testnoerror]("test 24", UPoPSObject:-BinaryAdd(upOne, UPoPSObject:-Negate(upOne)), 'assign'='upPmP');
Try[testnoerror]("test 25", UPoPSObject:-Truncate(upPmP), 'assign'='sub');
Try("test 26", evalb(expand(sub) = 0), true);

####################################################################

#Setup 
Try[testnoerror]("test 27",  UnivariatePolynomialOverPowerSeries(1), 'assign'='upOne');
Try[testnoerror]("test 28", UnivariatePolynomialOverPowerSeries(13*x^10*z + 2*x*y + 2*x + 1, 'x'), 'assign'='upP');
Try[testnoerror]("test 29", Truncate(upP), 'assign'='p');
Try[testnoerror]("test 30", UnivariatePolynomialOverPowerSeries(2*x*y + 13*z + 1, 'x'), 'assign'='upQ');
Try[testnoerror]("test 31", Truncate(upQ), 'assign'='q');

# Test Add
Try[testnoerror]("test 32",  Add(upP, upQ), 'assign'='upPaQ');
Try[testnoerror]("test 33", Truncate(upPaQ, 10), 'assign'='paq');
Try("test 34", evalb(expand(paq) = expand(p+q)), true);

Try[testnoerror]("test 35",  Add(upP, upP), 'assign'='upPaP');
Try[testnoerror]("test 36", Truncate(upPaP, 10), 'assign'='pap');
Try("test 37", evalb(expand(pap) = expand(2*p)), true);

Try[testnoerror]("test 38",  Add(upPaP, upPaQ), 'assign'='upPaPaPaQ');
Try[testnoerror]("test 39", Truncate(upPaPaPaQ, 10), 'assign'='papapaq');
Try("test 40", evalb(expand(papapaq) = expand(3*p + q)), true);

# Test Sub
Try[testnoerror]("test 41",  Subtract(upPaPaPaQ, upP), 'assign'='upPaPaQ');
Try[testnoerror]("test 42",  Subtract(upPaPaQ, upP), 'assign'='upPaQ');
Try[testnoerror]("test 43",  Subtract(upPaQ, upP), 'assign'='upaQ');
Try[testnoerror]("test 44",  Subtract(upaQ, upQ), 'assign'='upsub');
Try[testnoerror]("test 45", Truncate(upsub, 10), 'assign'='sub');
Try("test 46", evalb(expand(sub) = 0), true);

# Test Negate 
Try[testnoerror]("test 47",  Add(upP,  Negate(upP)), 'assign'='upPmP');
Try[testnoerror]("test 48", Truncate(upPmP), 'assign'='sub');
Try("test 49", evalb(expand(sub) = 0), true);

Try[testnoerror]("test 50",  Add(upOne,  Negate(upOne)), 'assign'='upPmP');
Try[testnoerror]("test 51", Truncate(upPmP), 'assign'='sub');
Try("test 52", evalb(expand(sub) = 0), true);

#end test