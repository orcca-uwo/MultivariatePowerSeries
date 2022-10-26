#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PuSO := MultivariatePowerSeries:-OldPuiseuxSeriesObject:
kernelopts(opaquemodules=true):

# Setup
pso := PowerSeries(2 + 3*(x + y) + 2*(x^2 + y^2) + x*y^2);
expDivs := table([x=2, y=3]);
Try[testnoerror]("test 1", PuSO(pso, expDivs), 'assign'='puso1');

pso := Inverse(pso);
expDivs :=table([x=2, y=1]);
Try[testnoerror]("test 2", PuSO(pso, expDivs), 'assign'='puso4');

pso := PowerSeries(2 + 2*(x + y));
expDivs := table([x=6, y=5]);
Try[testnoerror]("test 3", PuSO(pso, expDivs), 'assign'='puso2');

pso := Inverse(pso);
expDivs :=table([x=4, y=6]);
Try[testnoerror]("test 4", PuSO(pso, expDivs), 'assign'='puso3');

Try[testnoerror]("test 5", PuSO(PowerSeries(1), expDivs), 'assign'='pusoOne');

Try[testnoerror]("test 6", PuSO(PowerSeries(0), expDivs), 'assign'='pusoZero');

Try[testnoerror]("test 7", PuSO(PowerSeries(5), expDivs), 'assign'='pusoConst');


# Test Add
Try[testnoerror]("test 8", PuSO:-BinaryAdd(pusoZero, pusoZero), 'assign'='psAdd');

Try("test 9", psAdd:-ApproximatelyEqual(psAdd, pusoZero), true);

Try[testnoerror]("test 9", PuSO:-BinaryAdd(pusoOne, pusoZero), 'assign'='psAdd');

Try("test 10", psAdd:-ApproximatelyEqual(psAdd, pusoOne), true);

Try[testnoerror]("test 11", PuSO:-BinaryAdd(pusoOne, pusoConst), 'assign'='psAdd');

Try[testnoerror]("test 12", PuSO(PowerSeries(6), []), 'assign'='checkAdd');
Try("test 13", psAdd:-ApproximatelyEqual(psAdd, checkAdd), true);

Try[testnoerror]("test 14", PuSO:-BinaryAdd(puso1, pusoConst), 'assign'='psAdd');

pso := PowerSeries(7 + 3*(x + y) + 2*(x^2 + y^2) + x*y^2);
Try[testnoerror]("test 15", PuSO(pso, [x=2,y=3]), 'assign'='checkAdd');
Try("test 16", psAdd:-ApproximatelyEqual(psAdd, checkAdd, 3), true);

Try[testnoerror]("test 17", PuSO:-BinaryAdd(puso1, puso2), 'assign'='psAdd');

pso := PowerSeries(4 + 3*x^3 + 3*y^5 + 2*x^6 + 2*y^(10) + x^3*y^(10) + 2*x + 2*y^3);
Try[testnoerror]("test 18", PuSO(pso, [x=6,y=15]), 'assign'='checkAdd');
Try("test 19", psAdd:-ApproximatelyEqual(psAdd, checkAdd, 4), true);

Try[testnoerror]("test 20", PuSO:-BinaryAdd(puso3, puso4), 'assign'='psAdd');

pso := PowerSeries(1 - (1/2)*x - (1/2)*y^1 - (1/4)*x^2 + x*y + (1/2)*y^2
					+(5/8)*x^4 -(3/4)*y^6 + (9/4)*x^2*y^6 +(5/8)*y^12);
Try[testnoerror]("test 21", PuSO(pso, [x=4,y=6]), 'assign'='checkAdd');
Try("test 22", psAdd:-ApproximatelyEqual(psAdd, checkAdd, 2), true);

#################################################################################
# Setup
pso := PowerSeries(2 + 3*(x + y) + 2*(x^2 + y^2) + x*y^2 + z^2);
expDivs := table([x=2, y=3]);
Try[testnoerror]("test 23", PuSO(pso, expDivs), 'assign'='puso5');

pso := PowerSeries(2 + 2*(x + y)+ x^2*w);
expDivs := table([x=6, y=5, w=2]);
Try[testnoerror]("test 24", PuSO(pso, expDivs), 'assign'='puso6');

pso := PowerSeries(d -> x^d/d!, analytic=exp(x));
expDivs := table([x=2, y=5, w=4]);
Try[testnoerror]("test 25", PuSO(pso, expDivs), 'assign'='puso7');

# Test Add
Try[testnoerror]("test 26", PuSO:-BinaryAdd(puso5, puso6), 'assign'='psAdd');
Try[verify,table]("test 27", psAdd:-GetExpDivs(psAdd), 
								table(':-sparse' = 1, [x=6, y=15, z=1, w=2]));

pso := PowerSeries(4 + 3*x^3 + 3*y^5 + 2*x^6 + 2*y^(10) + x^3*y^(10) + z ^2+ 2*x + 2*y^3 + x^2*w);
Try[testnoerror]("test 28", PuSO(pso, [x=6,y=15, w=2, z=1]), 'assign'='checkAdd');
Try("test 29", psAdd:-ApproximatelyEqual(psAdd, checkAdd), true);								

Try[testnoerror]("test 30", PuSO:-BinaryAdd(puso6, puso7), 'assign'='psAdd');

pso := PowerSeries(3 + 2*x +2*y + x^2*w + x^3+ (1/2)*x^6 +(1/6)*x^9+(1/24)*x^12);
Try[testnoerror]("test 31", PuSO(pso, [x=6,y=5, w=2]), 'assign'='checkAdd');
Try("test 32", psAdd:-ApproximatelyEqual(psAdd, checkAdd, 12), true);	

Try[testnoerror]("test 33", PuSO:-BinaryAdd(puso3, puso7), 'assign'='psAdd');

Try[verify,table]("test 34", psAdd:-GetExpDivs(psAdd), 
								table(':-sparse' = 1, [x=4, y=6]));

pso := PowerSeries(3/2 -(1/2)*x -(1/2)*y + (3/2)*x^2 + x*y +(1/2)*y^2);
Try[testnoerror]("test 35", PuSO(pso, [x=4,y=6]), 'assign'='checkAdd');
Try("test 36", psAdd:-ApproximatelyEqual(psAdd, checkAdd, 2), true);	

pso := PowerSeries(d -> x^d/d!, analytic=exp(x));

Try[testnoerror]("test 37", PuSO:-BinaryAdd(puso7, pso), 'assign'='psAdd');

Try[verify,table]("test 38", psAdd:-GetExpDivs(psAdd), 
								table(':-sparse' = 1, [x=2]));

pso := PowerSeries(2 + x + (3/2)*x^2 + (1/6)*x^3 + (13/24)*x^4);
Try[testnoerror]("test 39", PuSO(pso, [x=2]), 'assign'='checkAdd');
Try("test 40", psAdd:-ApproximatelyEqual(psAdd, checkAdd, 4), true);	

pso := PowerSeries(d -> x^d/d!, analytic=exp(x));

Try[testnoerror]("test 41", PuSO:-BinaryAdd(pso, puso7), 'assign'='psAdd');

Try[verify,table]("test 42", psAdd:-GetExpDivs(psAdd), 
								table(':-sparse' = 1, [x=2]));

pso := PowerSeries(2 + x + (3/2)*x^2 + (1/6)*x^3 + (13/24)*x^4);
Try[testnoerror]("test 43", PuSO(pso, [x=2]), 'assign'='checkAdd');
Try("test 44", psAdd:-ApproximatelyEqual(psAdd, checkAdd, 4), true);	

#################################################################################
# Setup

upops := UnivariatePolynomialOverPowerSeries([PowerSeries(1), PowerSeries(0), PowerSeries(x), 
										  PowerSeries(y), 1/PowerSeries(1 + x + y)], z);

Try[testnoerror]("test 45", PuSO:-BinaryAdd(puso7, upops), 'assign'='psAdd');	

pso := PowerSeries(2 + x^(1) + (1/2)*x^2 + (1/6)*x^3 + x^2*z^2 + y*z^3 + z^4 +(1/24)*x^4);
Try[testnoerror]("test 46", PuSO(pso, [x=2, z=1]), 'assign'='checkAdd');
Try("test 47", psAdd:-ApproximatelyEqual(psAdd, checkAdd, 4), true);


Try[testnoerror]("test 48", PuSO:-BinaryAdd(upops, puso7), 'assign'='psAdd');	

Try("test 49", psAdd:-ApproximatelyEqual(psAdd, checkAdd, 4), true);
#end test