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
Try[testnoerror]("test 4", PSO:-Identity(x), 'assign'='psI');

# n-ary `+` operator 
Try[testnoerror]("test 4", psZero + psOne + psConst + psI, 'assign'='ps');
Try("test 5", PSO:-HomogeneousPart(ps, 0), 6);
Try("test 6", PSO:-HomogeneousPart(ps, 1), x);
Try("test 7", PSO:-HomogeneousPart(ps, 2), 0);
Try("test 8", PSO:-HomogeneousPart(ps, 10), 0);

Try[testnoerror]("test 9", 0 + psZero + psOne + 2 + 3*x + psConst  + x*y*z + psI + RootOf(x^2+1) , 'assign'='ps');
Try("test 10", PSO:-HomogeneousPart(ps, 0), 8 + RootOf(x^2+1));
Try("test 11", PSO:-HomogeneousPart(ps, 1), 4*x);
Try("test 12", PSO:-HomogeneousPart(ps, 3), x*y*z);
Try("test 13", PSO:-HomogeneousPart(ps, 10), 0);

Try[testnoerror]("test 14", PSO:-NaryAdd(psZero, psOne, psConst, psI, ':-coefficients' = [1, 2, 3, sqrt(2)]) , 'assign'='ps');
Try("test 15", PSO:-HomogeneousPart(ps, 0), 17);
Try[verify,simplify]("test 16", PSO:-HomogeneousPart(ps, 1), x*RootOf(_Z^2-2,index = 1));
Try("test 17", PSO:-HomogeneousPart(ps, 3), 0);
Try("test 18", PSO:-HomogeneousPart(ps, 10), 0);

Try[testnoerror]("test 19", PSO:-NaryAdd(psZero, psOne, psConst, psI, ':-coefficients' = [2, 2*x + y, 3*x*y, RootOf(x^2 + 1)]) , 'assign'='ps');
Try("test 20", PSO:-HomogeneousPart(ps, 0), 0);
Try("test 21", PSO:-HomogeneousPart(ps, 1), 2*x + y + RootOf(x^2+1)*x);
Try("test 22", PSO:-HomogeneousPart(ps, 2), 15*x*y);
Try("test 23", PSO:-HomogeneousPart(ps, 10), 0);


#end test
