#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PSO := MultivariatePowerSeries:-PowerSeriesObject:
kernelopts(opaquemodules=true):

# Setup 
Try[testnoerror]("test 1", PSO:-One(), 'assign'='psOne');
Try[testnoerror]("test 2", PSO:-Constant(5), 'assign'='psConst');
Try[testnoerror]("test 3", PSO:-Identity(x), 'assign'='psI');

# n-ary `*` operator 
Try[testnoerror]("test 4", PSO:-`*`(psOne, psConst, psI, 1, -1, -1), 'assign'='ps');
Try("test 5", evalb(PSO:-Truncate(ps, 5) = 5*x), true);
Try("test 6", PSO:-HomogeneousPart(ps, 1), 5*x);
Try("test 7", PSO:-ApproximatelyZero(ps, 'force'), false);

Try[testnoerror]("test 8", PSO:-GeometricSeries([x, y, z]), 'assign'='ps1');
Try[testnoerror]("test 9", PSO:-SumOfAllMonomials([x, y, z]), 'assign'='ps2');
Try("test 10", PSO:-HomogeneousPart(PSO:-`*`(ps1, ps2), 0), 1);
Try("test 11", PSO:-HomogeneousPart(PSO:-`*`(ps1, ps2), 1), 2*x + 2*y + 2*z);
Try("test 12", degree(HomogeneousPart(PSO:-`*`(ps1, ps2), 8)), 8);

#end test
