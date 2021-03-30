#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PSO := MultivariatePowerSeries:-PowerSeriesObject:
kernelopts(opaquemodules=true):

# Setup
Try[testnoerror]("test 1", PSO:-Zero(), 'assign'='psZero');
Try[testnoerror]("test 2", PSO:-One(), 'assign'='psOne');
Try[testnoerror]("test 3", PSO:-Constant(5), 'assign'='psConst');
Try[testnoerror]("test 4", PSO:-FromPolynomial(2 + 1/3*(x + y)), 'assign'='psP');
Try[testnoerror]("test 5", PSO:-Truncate(psP), 'assign'='p');
Try[testnoerror]("test 6", PSO:-FromPolynomial(2 + 1/3*(x + y) + 2/3*(x^2 + y^2) + x*y^2), 'assign'='psQ');
Try[testnoerror]("test 7", PSO:-Truncate(psQ), 'assign'='q');

# Test Add
Try[testnoerror]("test 8", PSO:-BinaryAdd(psZero, psZero), 'assign'='psAdd');
Try("test 9", PSO:-ApproximatelyZero(psAdd), true);

Try[testnoerror]("test 10", PSO:-BinaryAdd(psOne, psZero), 'assign'='psAdd');
Try("test 11", PSO:-ApproximatelyZero(psAdd), false);
Try("test 12", PSO:-HomogeneousPart(psAdd, 0), 1);
Try("test 13", PSO:-HomogeneousPart(psAdd, 1), 0);
Try("test 14", PSO:-HomogeneousPart(psAdd, 10), 0);

Try[testnoerror]("test 15", PSO:-BinaryAdd(psZero, psOne), 'assign'='psAdd');
Try("test 16", PSO:-ApproximatelyZero(psAdd), false);
Try("test 17", PSO:-HomogeneousPart(psAdd, 0), 1);
Try("test 18", PSO:-HomogeneousPart(psAdd, 1), 0);
Try("test 19", PSO:-HomogeneousPart(psAdd, 10), 0);

Try[testnoerror]("test 20", PSO:-BinaryAdd(psConst, psOne), 'assign'='psAdd');
Try("test 21", PSO:-HomogeneousPart(psAdd, 0), 6);
Try("test 22", PSO:-HomogeneousPart(psAdd, 1), 0);

Try[testnoerror]("test 23", PSO:-BinaryAdd(psP, psP), 'assign'='psAdd');
Try[testnoerror]("test 24", PSO:-Truncate(psAdd), 'assign'='pp');
Try("test 25", evalb(pp = expand(2*p)), true);

Try[testnoerror]("test 26", PSO:-BinaryAdd(psAdd, psQ), 'assign'='psAdd');
Try[testnoerror]("test 27", PSO:-Truncate(psAdd, 10), 'assign'='ppq');
Try("test 28", evalb(ppq = expand(pp + q)), true);

# Test Algebraic
Try[testnoerror]("test 29", PSO:-BinaryAdd(PSO:-Constant(RootOf(z^2 + 1)), PSO:-Constant(RootOf(z^2 + 1))), 'assign'='psAddAlg');
Try("test 30", PSO:-ApproximatelyZero(psAddAlg), false);
Try("test 31", PSO:-HomogeneousPart(psAddAlg, 0), 2*RootOf(z^2+1) );
Try("test 32", PSO:-HomogeneousPart(psAddAlg, 1), 0);


# Test Sub
Try[testnoerror]("test 33", PSO:-BinarySub(psAdd, psQ), 'assign'='psSub');
Try[testnoerror]("test 34", PSO:-Truncate(psSub), 'assign'='pp');
Try("test 35", evalb(pp = expand(2*p)), true);

Try[testnoerror]("test 36", PSO:-BinarySub(psSub, psP), 'assign'='psSub');
Try[testnoerror]("test 37", PSO:-Truncate(psSub), 'assign'='pps');
Try("test 38", evalb(pps = p), true);

Try[testnoerror]("test 39", PSO:-BinarySub(psSub, psP), 'assign'='psSub');
Try[testnoerror]("test 40", PSO:-Truncate(psSub), 'assign'='pss');
Try("test 41", evalb(pss = 0), true);

# Test Negate 
Try[testnoerror]("test 42", PSO:-Negate(psOne), 'assign'='psNOne');
Try[testnoerror]("test 43", PSO:-Truncate(psNOne), 'assign'='neg');
Try("test 44", evalb(neg = -1), true);

Try[testnoerror]("test 45", PSO:-BinaryAdd(psP, PSO:-Negate(psP)), 'assign'='psNP');
Try[testnoerror]("test 46", PSO:-Truncate(psNP), 'assign'='neg');
Try("test 47", evalb(neg = 0), true);

Try[testnoerror]("test 48", PSO:-BinaryAdd(PSO:-Negate(psQ), psQ), 'assign'='psNQ');
Try[testnoerror]("test 49", PSO:-Truncate(psNQ), 'assign'='neg');
Try("test 50", evalb(neg = 0), true);

# Test Algebraic
Try[testnoerror]("test 51", PSO:-BinaryAdd(PSO:-Constant(RootOf(z^2 + 1)), PSO:-Negate(PSO:-Constant(RootOf(z^2 + 1)))), 'assign'='psAddAlg');
Try("test 52", PSO:-ApproximatelyZero(psAddAlg), true);
Try("test 53", PSO:-HomogeneousPart(psAddAlg, 0), 0);
Try("test 54", PSO:-HomogeneousPart(psAddAlg, 1), 0);


####################################################################
# Setup
Try[testnoerror]("test 55", PowerSeries(0), 'assign'='psZero');
Try[testnoerror]("test 56", PowerSeries(1), 'assign'='psOne');
Try[testnoerror]("test 57", PowerSeries(5), 'assign'='psConst');
Try[testnoerror]("test 58", PowerSeries(2 + 1/3*(x + y)), 'assign'='psP');
Try[testnoerror]("test 59", Truncate(psP), 'assign'='p');
Try[testnoerror]("test 60", PowerSeries(2 + 1/3*(x + y) + 2/3*(x^2 + y^2) + x*y^2), 'assign'='psQ');
Try[testnoerror]("test 61", Truncate(psQ), 'assign'='q');

# Test Add
Try[testnoerror]("test 62",  Add(psZero, psZero), 'assign'='psAdd');
Try("test 63",  ApproximatelyZero(psAdd), true);

Try[testnoerror]("test 64",  Add(psOne, psZero), 'assign'='psAdd');
Try("test 65",  ApproximatelyZero(psAdd), false);
Try("test 66",  HomogeneousPart(psAdd, 0), 1);
Try("test 67",  HomogeneousPart(psAdd, 1), 0);
Try("test 68",  HomogeneousPart(psAdd, 10), 0);

Try[testnoerror]("test 69",  Add(psZero, psOne), 'assign'='psAdd');
Try("test 70",  ApproximatelyZero(psAdd), false);
Try("test 71",  HomogeneousPart(psAdd, 0), 1);
Try("test 72",  HomogeneousPart(psAdd, 1), 0);
Try("test 73",  HomogeneousPart(psAdd, 10), 0);

Try[testnoerror]("test 74",  Add(psConst, psOne), 'assign'='psAdd');
Try("test 75",  HomogeneousPart(psAdd, 0), 6);
Try("test 76",  HomogeneousPart(psAdd, 1), 0);

Try[testnoerror]("test 77",  Add(psP, psP), 'assign'='psAdd');
Try[testnoerror]("test 78", Truncate(psAdd), 'assign'='pp');
Try("test 79", evalb(pp = expand(2*p)), true);

Try[testnoerror]("test 80",  Add(psAdd, psQ), 'assign'='psAdd');
Try[testnoerror]("test 81", Truncate(psAdd, 10), 'assign'='ppq');
Try("test 82", evalb(ppq = expand(pp + q)), true);

# Test Algebraic
Try[testnoerror]("test 83",  Add(PowerSeries(RootOf(z^2 + 1)), PowerSeries(RootOf(z^2 + 1))), 'assign'='psAddAlg');
Try("test 84",  ApproximatelyZero(psAddAlg), false);
Try("test 85",  HomogeneousPart(psAddAlg, 0), 2*RootOf(z^2+1) );
Try("test 86",  HomogeneousPart(psAddAlg, 1), 0);


# Test Sub
Try[testnoerror]("test 87",  Subtract(psAdd, psQ), 'assign'='psSub');
Try[testnoerror]("test 88", Truncate(psSub), 'assign'='pp');
Try("test 89", evalb(pp = expand(2*p)), true);

Try[testnoerror]("test 90",  Subtract(psSub, psP), 'assign'='psSub');
Try[testnoerror]("test 91", Truncate(psSub), 'assign'='pps');
Try("test 92", evalb(pps = p), true);

Try[testnoerror]("test 93",  Subtract(psSub, psP), 'assign'='psSub');
Try[testnoerror]("test 94", Truncate(psSub), 'assign'='pss');
Try("test 95", evalb(pss = 0), true);

# Test Negate 
Try[testnoerror]("test 96",  Negate(psOne), 'assign'='psNOne');
Try[testnoerror]("test 97", Truncate(psNOne), 'assign'='neg');
Try("test 98", evalb(neg = -1), true);

Try[testnoerror]("test 99",  Add(psP,  Negate(psP)), 'assign'='psNP');
Try[testnoerror]("test 100", Truncate(psNP), 'assign'='neg');
Try("test 101", evalb(neg = 0), true);

Try[testnoerror]("test 102",  Add( Negate(psQ), psQ), 'assign'='psNQ');
Try[testnoerror]("test 103", Truncate(psNQ), 'assign'='neg');
Try("test 104", evalb(neg = 0), true);

# Test Algebraic
Try[testnoerror]("test 105",  Add(PowerSeries(RootOf(z^2 + 1)),  Negate(PowerSeries(RootOf(z^2 + 1)))), 'assign'='psAddAlg');
Try("test 106",  ApproximatelyZero(psAddAlg), true);
Try("test 107",  HomogeneousPart(psAddAlg, 0), 0);
Try("test 108",  HomogeneousPart(psAddAlg, 1), 0);

#end test
