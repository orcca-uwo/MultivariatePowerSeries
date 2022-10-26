#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));
kernelopts(opaquemodules=false):
UPoPSObject := MultivariatePowerSeries:-UnivariatePolynomialOverPowerSeriesObject:
kernelopts(opaquemodules=true):

#Setup
Try[testnoerror]("test 1", UPoPSObject:-Zero(), 'assign'='upZero');

# Test TaylorShift
Try[testnoerror]("test 2", UPoPSObject:-TaylorShift(upZero, 0), 'assign'='upTS');
Try("test 3", UPoPSObject:-ApproximatelyEqual(upZero, upTS, ':-force'), true);

Try[testnoerror]("test 4", UPoPSObject:-TaylorShift(upZero, 1), 'assign'='upTS');
Try("test 5", UPoPSObject:-ApproximatelyEqual(upZero, upTS, ':-force'), true);

Try[testnoerror]("test 6", UPoPSObject:-TaylorShift(upZero, 0), 'assign'='upTS');
Try("test 7", UPoPSObject:-ApproximatelyEqual(upZero, upTS, ':-force'), true);

Try[testnoerror]("test 8", UPoPSObject:-FromPolynomial(13*x^10*z - 3*y*x^5 + 2*x*y + 2*x + 1, 'x'), 'assign'='upP');
Try[testnoerror]("test 9", UPoPSObject:-TaylorShift(upP, 0), 'assign'='upTS');
Try[testnoerror]("test 10", UPoPSObject:-Truncate(upTS, 10), 'assign'='ts');
Try("test 11", UPoPSObject:-ApproximatelyEqual(upTS, upP, ':-force'), true);
Try[testnoerror]("test 12", x -> expand(13*x^10*z - 3*y*x^5 + 2*x*y + 2*x + 1) , 'assign'='func');
Try("test 13", evalb(expand(ts) = func(x + 0)), true);

Try[testnoerror]("test 14", UPoPSObject:-TaylorShift(upP, 10), 'assign'='upTS');
Try[testnoerror]("test 15", UPoPSObject:-Truncate(upTS, 10), 'assign'='ts');
Try("test 16", evalb(expand(ts) = func(x + 10)), true);

Try[testnoerror]("test 17", UPoPSObject:-TaylorShift(upP, -10/3), 'assign'='upTS');
Try[testnoerror]("test 18", UPoPSObject:-Truncate(upTS, 10), 'assign'='ts');
Try("test 19", evalb(expand(ts) = func(x - 10/3)), true);

Try[testnoerror]("test 20", randpoly([x, y, z], terms = 10), 'assign'='poly');
Try[testnoerror]("test 21", X -> eval(poly, x = X), 'assign'='func');
Try[testnoerror]("test 22", UPoPSObject:-FromPolynomial(poly, 'x'), 'assign'='upP');
Try[testnoerror]("test 23",  UPoPSObject:-TaylorShift(upP, RootOf(x^2+1)), 'assign'='upTS');
Try[testnoerror]("test 24", UPoPSObject:-Truncate(upTS, 20), 'assign'='ts');
Try("test 25", evalb(expand(ts) = Algebraic:-Expand(func(x + RootOf(x^2+1)))), true);

####################################################################

#Setup
Try[testnoerror]("test 26", UnivariatePolynomialOverPowerSeries(0), 'assign'='upZero');

# Test TaylorShift
Try[testnoerror]("test 27", TaylorShift(upZero, 0), 'assign'='upTS');
Try("test 28", ApproximatelyEqual(upZero, upTS, ':-force'), true);

Try[testnoerror]("test 29", TaylorShift(upZero, 1), 'assign'='upTS');
Try("test 30", ApproximatelyEqual(upZero, upTS, ':-force'), true);

Try[testnoerror]("test 31", TaylorShift(upZero, 0), 'assign'='upTS');
Try("test 32", ApproximatelyEqual(upZero, upTS, ':-force'), true);

Try[testnoerror]("test 33", UnivariatePolynomialOverPowerSeries(13*x^10*z - 3*y*x^5 + 2*x*y + 2*x + 1, 'x'), 'assign'='upP');
Try[testnoerror]("test 34", TaylorShift(upP, 0), 'assign'='upTS');
Try[testnoerror]("test 35", Truncate(upTS, 10), 'assign'='ts');
Try("test 36", ApproximatelyEqual(upTS, upP, ':-force'), true);
Try[testnoerror]("test 37", x -> expand(13*x^10*z - 3*y*x^5 + 2*x*y + 2*x + 1) , 'assign'='func');
Try("test 38", evalb(expand(ts) = func(x + 0)), true);

Try[testnoerror]("test 39", TaylorShift(upP, 10), 'assign'='upTS');
Try[testnoerror]("test 40", Truncate(upTS, 10), 'assign'='ts');
Try("test 41", evalb(expand(ts) = func(x + 10)), true);

Try[testnoerror]("test 42", TaylorShift(upP, -10/3), 'assign'='upTS');
Try[testnoerror]("test 43", Truncate(upTS, 10), 'assign'='ts');
Try("test 44", evalb(expand(ts) = func(x - 10/3)), true);


Try[testnoerror]("test 45", randpoly([x, y, z], terms = 10), 'assign'='poly');
Try[testnoerror]("test 46", X -> eval(poly, x = X), 'assign'='func');
Try[testnoerror]("test 47", UnivariatePolynomialOverPowerSeries(poly, 'x'), 'assign'='upP');
Try[testnoerror]("test 48",  TaylorShift(upP, RootOf(x^2+1)), 'assign'='upTS');
Try[testnoerror]("test 49", Truncate(upTS, 20), 'assign'='ts');
Try("test 50", evalb(expand(ts) = Algebraic:-Expand(func(x + RootOf(x^2+1)))), true);


#end test
