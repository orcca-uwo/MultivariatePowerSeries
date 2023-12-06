#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PSO := MultivariatePowerSeries:-PowerSeriesObject:
kernelopts(opaquemodules=true):

# Setup
Try[testnoerror]("test 1", PSO:-Zero(), 'assign'='psZero');
Try[testnoerror]("test 2", PSO:-One(), 'assign'='psOne');
Try[testnoerror]("test 3", PSO:-FromPolynomial(2 + 1/3*(x + y)), 'assign'='psP');
Try[testnoerror]("test 4", PSO:-Truncate(psP), 'assign'='p');
Try[testnoerror]("test 5", PSO:-FromPolynomial(2 + 1/3*(x + y) + 2/3*(x^2 + y^2) + x*y^2), 'assign'='psQ');
Try[testnoerror]("test 6", PSO:-Truncate(psQ), 'assign'='q');

# Test Product
Try[testnoerror]("test 7", PSO:-BinaryMultiply(psOne, psZero), 'assign'='psProd');
Try("test 8", PSO:-ApproximatelyZero(psProd, 'force'), true);
Try("test 9", PSO:-HomogeneousPart(psProd, 1), 0);

Try[testnoerror]("test 10", PSO:-BinaryMultiply(psZero, psP), 'assign'='psProd');
Try("test 11", PSO:-ApproximatelyZero(psProd, 'force'), true);
Try("test 12", PSO:-HomogeneousPart(psProd, 1), 0);

Try[testnoerror]("test 13", PSO:-BinaryMultiply(psQ, psP), 'assign'='psProd');
Try[testnoerror]("test 14", PSO:-Truncate(psProd, 10), 'assign'='pq');
Try("test 15", evalb(pq = expand(p*q)), true);

Try[testnoerror]("test 16", PSO:-BinaryMultiply(psOne, psP), 'assign'='psProd');
Try[testnoerror]("test 17", PSO:-Truncate(psProd, 5), 'assign'='ppp');
Try("test 18", evalb(ppp = p), true);

# Test Algebraic
Try[testnoerror]("test 19", PSO:-FromPolynomial(RootOf(z^2 + 1)^2 * x + RootOf(z^2+1)), 'assign'='psAlg1');
Try[testnoerror]("test 20", PSO:-Truncate(psAlg1), 'assign'='alg1');
Try("test 21", evalb( alg1 = RootOf(z^2+1) - x ), true);


Try[testnoerror]("test 22", PSO:-FromPolynomial(RootOf(z^2 + 1)*sqrt(2)*x^2 + sqrt(3)*x + x*y + RootOf(z^2 + 1)), 'assign'='psAlg2');
Try[testnoerror]("test 23", PSO:-Truncate(psAlg2), 'assign'='alg2');
Try("test 24", evalb( alg2 = RootOf(_Z^2 + 1)*sqrt(2)*x^2 + sqrt(3)*x + x*y + RootOf(_Z^2 + 1) ), true);

Try[testnoerror]("test 25", PSO:-BinaryMultiply(psAlg1, psAlg2), 'assign'='psProd');
Try("test 26", PSO:-ApproximatelyZero(psProd, 'force'), false);
Try("test 27", PSO:-HomogeneousPart(psProd, 0), -1);
Try[verify,simplify]("test 28", PSO:-HomogeneousPart(psProd, 1), RootOf(_Z^2+1)*RootOf(_Z^2-3,index = 1)*x-x*RootOf(_Z^2+1));
Try[verify,simplify]("test 29", PSO:-HomogeneousPart(psProd, 2), -x^2*RootOf(_Z^2-2,index = 1)-x^2*RootOf(_Z^2-3,index = 1)+RootOf(_Z^2+1)*x*y);
Try[verify,simplify]("test 30", PSO:-HomogeneousPart(psProd, 3), -RootOf(_Z^2+1)*RootOf(_Z^2-2,index = 1)*x^3-x^2*y);
Try("test 31", PSO:-HomogeneousPart(psProd, 4), 0);


####################################################################

# Setup
Try[testnoerror]("test 0.0", PowerSeries(0), 'assign'='psZero');
Try[testnoerror]("test 1.0", PowerSeries(1), 'assign'='psOne');
Try[testnoerror]("test 4.17", PowerSeries(2 + 1/3*(x + y)), 'assign'='psP');
Try[testnoerror]("test 4.28", Truncate(psP), 'assign'='p');
Try[testnoerror]("test 4.30", PowerSeries(2 + 1/3*(x + y) + 2/3*(x^2 + y^2) + x*y^2), 'assign'='psQ');
Try[testnoerror]("test 4.41", Truncate(psQ), 'assign'='q');

# Test Product
Try[testnoerror]("test 11.0", Multiply(psOne, psZero), 'assign'='psProd');
Try("test 11.1", ApproximatelyZero(psProd, 'force'), true);
Try("test 11.2", HomogeneousPart(psProd, 1), 0);

Try[testnoerror]("test 11.3", Multiply(psZero, psP), 'assign'='psProd');
Try("test 11.4", ApproximatelyZero(psProd, 'force'), true);
Try("test 11.5", HomogeneousPart(psProd, 1), 0);

Try[testnoerror]("test 11.6", Multiply(psQ, psP), 'assign'='psProd');
Try[testnoerror]("test 11.7", Truncate(psProd, 10), 'assign'='pq');
Try("test 11.8", evalb(pq = expand(p*q)), true);

Try[testnoerror]("test 11.9", Multiply(psOne, psP), 'assign'='psProd');
Try[testnoerror]("test 11.10", Truncate(psProd, 5), 'assign'='ppp');
Try("test 11.11", evalb(ppp = p), true);

# Test Algebraic
Try[testnoerror]("test 11.11", PowerSeries(RootOf(z^2 + 1)^2 * x + RootOf(z^2+1)), 'assign'='psAlg1');
Try[testnoerror]("test 11.12", Truncate(psAlg1), 'assign'='alg1');
Try("test 11.13", evalb( alg1 = RootOf(z^2+1) - x ), true);


Try[testnoerror]("test 11.14", PowerSeries(RootOf(z^2 + 1)*sqrt(2)*x^2 + sqrt(3)*x + x*y + RootOf(z^2 + 1)), 'assign'='psAlg2');
Try[testnoerror]("test 11.15", Truncate(psAlg2), 'assign'='alg2');
Try("test 11.16", evalb( alg2 = RootOf(_Z^2 + 1)*sqrt(2)*x^2 + sqrt(3)*x + x*y + RootOf(_Z^2 + 1) ), true);

Try[testnoerror]("test 11.17", BinaryMultiply(psAlg1, psAlg2), 'assign'='psProd');
Try("test 11.18", ApproximatelyZero(psProd, 'force'), false);
Try("test 11.20", HomogeneousPart(psProd, 0), -1);
Try[verify,simplify]("test 11.21", HomogeneousPart(psProd, 1), RootOf(_Z^2+1)*RootOf(_Z^2-3,index = 1)*x-x*RootOf(_Z^2+1));
Try[verify,simplify]("test 11.22", HomogeneousPart(psProd, 2), -x^2*RootOf(_Z^2-2,index = 1)-x^2*RootOf(_Z^2-3,index = 1)+RootOf(_Z^2+1)*x*y);
Try[verify,simplify]("test 11.23", HomogeneousPart(psProd, 3), -RootOf(_Z^2+1)*RootOf(_Z^2-2,index = 1)*x^3-x^2*y);
Try("test 11.24", HomogeneousPart(psProd, 4), 0);


#end test