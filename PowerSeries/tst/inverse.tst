#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PSO := MultivariatePowerSeries:-PowerSeriesObject:
kernelopts(opaquemodules=true):

# Setup
Try[testnoerror]("test 1", PSO:-Constant(5), 'assign'='psConst');
Try[testnoerror]("test 2", PSO:-FromPolynomial(2 + 1/3*(x + y)), 'assign'='psP');
Try[testnoerror]("test 3", PSO:-FromPolynomial(2 + 1/3*(x + y) + 2/3*(x^2 + y^2) + x*y^2), 'assign'='psQ');

# Test ExactQuotient and Inverse
Try[testnoerror]("test 4", PSO:-Inverse(psP), 'assign'='psPInv');
Try[testnoerror]("test 5", PSO:-BinaryMultiply(psP, psPInv), 'assign'='psPPinv');
Try[testnoerror]("test 6", PSO:-Truncate(psPPinv, 10), 'assign'='ppinv');
Try("test 7", evalb(ppinv = 1), true);

Try[testnoerror]("test 8", PSO:-Inverse(psQ), 'assign'='psQInv');
Try[testnoerror]("test 9", PSO:-BinaryMultiply(psQ, psQInv), 'assign'='psQQinv');
Try[testnoerror]("test 10", PSO:-Truncate(psQQinv, 15), 'assign'='qqinv');
Try("test 11", evalb(qqinv = 1), true);

Try[testnoerror]("test 12", PSO:-Inverse(psConst), 'assign'='psCInv');
Try[testnoerror]("test 13", PSO:-BinaryMultiply(psConst, psCInv), 'assign'='psCCinv');
Try[testnoerror]("test 14", PSO:-Truncate(psCCinv, 5), 'assign'='ccinv');
Try("test 15", evalb(ccinv = 1), true);

# Test Algebraic
Try[testnoerror]("test 16", PSO:-FromPolynomial(RootOf(z^2 + 1)*sqrt(2)*x^2 + x*y + sqrt(3)*x + sqrt(2)), 'assign'='psAlg2');
Try[testnoerror]("test 17", PSO:-Truncate(psAlg2), 'assign'='alg2');
Try("test 18", evalb( alg2 = RootOf(z^2 + 1)*sqrt(2)*x^2 + x*y + sqrt(3)*x + sqrt(2)), true);
Try[testnoerror]("test 19", PSO:-Inverse(psAlg2), 'assign'='psAlgInv');
Try[testnoerror]("test 20", PSO:-BinaryMultiply(psAlgInv, psAlg2), 'assign'='psAlgOne');
Try[testnoerror]("test 21", PSO:-Truncate(psAlgOne, 10), 'assign'='psao');
Try("test 22", evalb(psao = 1), true);


####################################################################

# Setup
Try[testnoerror]("test 23", PowerSeries(5), 'assign'='psConst');
Try[testnoerror]("test 24", PowerSeries(2 + 1/3*(x + y)), 'assign'='psP');
Try[testnoerror]("test 25", PowerSeries(2 + 1/3*(x + y) + 2/3*(x^2 + y^2) + x*y^2), 'assign'='psQ');

# Test ExactQuotient and Inverse
Try[testnoerror]("test 26", Inverse(psP), 'assign'='psPInv');
Try[testnoerror]("test 27", Multiply(psP, psPInv), 'assign'='psPPinv');
Try[testnoerror]("test 28", Truncate(psPPinv, 10), 'assign'='ppinv');
Try("test 29", evalb(ppinv = 1), true);

Try[testnoerror]("test 30", Inverse(psQ), 'assign'='psQInv');
Try[testnoerror]("test 31", Multiply(psQ, psQInv), 'assign'='psQQinv');
Try[testnoerror]("test 32", Truncate(psQQinv, 15), 'assign'='qqinv');
Try("test 33", evalb(qqinv = 1), true);

Try[testnoerror]("test 34", Inverse(psConst), 'assign'='psCInv');
Try[testnoerror]("test 35", Multiply(psConst, psCInv), 'assign'='psCCinv');
Try[testnoerror]("test 36", Truncate(psCCinv, 5), 'assign'='ccinv');
Try("test 37", evalb(ccinv = 1), true);

# Test Algebraic
Try[testnoerror]("test 38", PowerSeries(RootOf(z^2 + 1)*sqrt(2)*x^2 + x*y + sqrt(3)*x + sqrt(2)), 'assign'='psAlg2');
Try[testnoerror]("test 39", Truncate(psAlg2), 'assign'='alg2');
Try("test 40", evalb( alg2 = RootOf(z^2 + 1)*sqrt(2)*x^2 + x*y + sqrt(3)*x + sqrt(2) ), true);
Try[testnoerror]("test 41", Inverse(psAlg2), 'assign'='psAlgInv');
Try[testnoerror]("test 42", Multiply(psAlgInv, psAlg2), 'assign'='psAlgOne');
Try[testnoerror]("test 43", Truncate(psAlgOne, 10), 'assign'='psao');
Try("test 44", evalb(psao = 1), true);

#end test