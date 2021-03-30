#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PSO := MultivariatePowerSeries:-PowerSeriesObject:
kernelopts(opaquemodules=true):

# Test (From/To)Polynomial
Try[testnoerror]("test 0", PSO:-FromPolynomial(0), 'assign'='psP');
Try("test 1", PSO:-ApproximatelyZero(psP), true);
Try("test 2", PSO:-ApproximatelyEqual(psP, PSO:-One()), false);
Try("test 4", PSO:-ApproximatelyEqual(psP, psP), true);
Try("test 5", PSO:-HomogeneousPart(psP, 0), 0);
Try("test 6", PSO:-HomogeneousPart(psP, 10), 0);
Try[testnoerror]("test 7", PSO:-Truncate(psP), 'assign'='p');
Try("test 8", evalb(p = 0), true);

Try[testnoerror]("test 9", PSO:-FromPolynomial(1), 'assign'='psP');
Try("test 10", PSO:-ApproximatelyZero(psP), false);
Try("test 11", PSO:-ApproximatelyEqual(psP, PSO:-One()), true);
Try("test 12", PSO:-ApproximatelyEqual(psP, psP), true);
Try("test 13", PSO:-HomogeneousPart(psP, 0), 1);
Try("test 14", PSO:-HomogeneousPart(psP, 10), 0);
Try[testnoerror]("test 4.15", PSO:-Truncate(psP), 'assign'='p');
Try("test 16", evalb(p = 1), true);

Try[testnoerror]("test 4.17", PSO:-FromPolynomial(2 + 1/3*(x + y)), 'assign'='psP');
Try("test 18", PSO:-ApproximatelyZero(psP), false);
Try("test 19", PSO:-IsUnit(psP), true);
Try("test 20", PSO:-ApproximatelyEqual(psP, PSO:-One()), false);
Try("test 21", PSO:-ApproximatelyEqual(psP, psP), true);
Try("test 22", PSO:-HomogeneousPart(psP, 0), 2);
Try("test 23", PSO:-HomogeneousPart(psP, 1), x/3 + y/3);
Try("test 24", PSO:-HomogeneousPart(psP, 2), 0);
Try("test 25", PSO:-GetCoefficient(psP, 1), 2);
Try("test 26", PSO:-GetCoefficient(psP, x), 1/3);
Try("test 27", PSO:-GetCoefficient(psP, x^2*y^2), 0);
Try[testnoerror]("test 28", PSO:-Truncate(psP), 'assign'='p');
Try("test 29", evalb(p = 2 + 1/3*(x + y)), true);

Try[testnoerror]("test 30", PSO:-FromPolynomial(2 + 1/3*(x + y) + 2/3*(x^2 + y^2) + x*y^2), 'assign'='psQ');
Try("test 31", PSO:-ApproximatelyZero(psQ), false);
Try("test 32", PSO:-IsUnit(psQ), true);
Try("test 33", PSO:-ApproximatelyEqual(psQ, PSO:-One()), false);
Try("test 34", PSO:-ApproximatelyEqual(psQ, psP, 1), true);
Try("test 35", PSO:-ApproximatelyEqual(psQ, psP), false);
Try("test 36", PSO:-ApproximatelyEqual(psQ, psQ), true);
Try("test 37", PSO:-HomogeneousPart(psQ, 2), (2*x^2)/3 + (2*y^2)/3);
Try("test 38", PSO:-HomogeneousPart(psQ, 3), x*y^2);
Try("test 39", PSO:-HomogeneousPart(psQ, 4), 0);
Try("test 40", PSO:-GetCoefficient(psQ, x^2), 2/3);
Try[testnoerror]("test 41", PSO:-Truncate(psQ), 'assign'='q');
Try("test 42", evalb(q = 2 + 1/3*(x + y) + 2/3*(x^2 + y^2) + x*y^2), true);


####################################################################

Try[testnoerror]("test 43", PowerSeries(0), 'assign'='psP');
Try("test 44", ApproximatelyZero(psP), true);
Try("test 45", ApproximatelyEqual(psP, PowerSeries(1)), false);
Try("test 46", ApproximatelyEqual(psP, psP), true);
Try("test 47", HomogeneousPart(psP, 0), 0);
Try("test 48", HomogeneousPart(psP, 10), 0);
Try[testnoerror]("test 49", Truncate(psP), 'assign'='p');
Try("test 50", evalb(p = 0), true);

Try[testnoerror]("test 51", PowerSeries(1), 'assign'='psP');
Try("test 52", ApproximatelyZero(psP), false);
Try("test 53", ApproximatelyEqual(psP, PowerSeries(1)), true);
Try("test 54", ApproximatelyEqual(psP, psP), true);
Try("test 55", HomogeneousPart(psP, 0), 1);
Try("test 56", HomogeneousPart(psP, 10), 0);
Try[testnoerror]("test 57", Truncate(psP), 'assign'='p');
Try("test 58", evalb(p = 1), true);

Try[testnoerror]("test 59", PowerSeries(2 + 1/3*(x + y)), 'assign'='psP');
Try("test 60", ApproximatelyZero(psP), false);
Try("test 61", IsUnit(psP), true);
Try("test 62", ApproximatelyEqual(psP, PowerSeries(1)), false);
Try("test 63", ApproximatelyEqual(psP, psP), true);
Try("test 64", HomogeneousPart(psP, 0), 2);
Try("test 65", HomogeneousPart(psP, 1), x/3 + y/3);
Try("test 66", HomogeneousPart(psP, 2), 0);
Try("test 67", GetCoefficient(psP, 1), 2);
Try("test 68", GetCoefficient(psP, x), 1/3);
Try("test 69", GetCoefficient(psP, x^2*y^2), 0);
Try[testnoerror]("test 70", Truncate(psP), 'assign'='p');
Try("test 71", evalb(p = 2 + 1/3*(x + y)), true);

Try[testnoerror]("test 72", PowerSeries(2 + 1/3*(x + y) + 2/3*(x^2 + y^2) + x*y^2), 'assign'='psQ');
Try("test 73", ApproximatelyZero(psQ), false);
Try("test 74", IsUnit(psQ), true);
Try("test 75", ApproximatelyEqual(psQ, PowerSeries(1)), false);
Try("test 76", ApproximatelyEqual(psQ, psP, 1), true);
Try("test 77", ApproximatelyEqual(psQ, psP), false);
Try("test 78", ApproximatelyEqual(psQ, psQ), true);
Try("test 79", HomogeneousPart(psQ, 2), (2*x^2)/3 + (2*y^2)/3);
Try("test 80", HomogeneousPart(psQ, 3), x*y^2);
Try("test 81", HomogeneousPart(psQ, 4), 0);
Try("test 82", GetCoefficient(psQ, x^2), 2/3);
Try[testnoerror]("test 83", Truncate(psQ), 'assign'='q');
Try("test 84", evalb(q = 2 + 1/3*(x + y) + 2/3*(x^2 + y^2) + x*y^2), true);

#end test