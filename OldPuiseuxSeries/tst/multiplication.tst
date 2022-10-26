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
Try[testnoerror]("test 8", PuSO:-BinaryMultiply(pusoZero, pusoZero), 'assign'='psMult');

Try("test 9", psMult:-ApproximatelyEqual(psMult, pusoZero), true);

Try[testnoerror]("test 9", PuSO:-BinaryMultiply(pusoOne, pusoZero), 'assign'='psMult');

Try("test 10", psMult:-ApproximatelyEqual(psMult, pusoZero), true);

Try[testnoerror]("test 11", PuSO:-BinaryMultiply(pusoOne, pusoConst), 'assign'='psMult');

Try("test 13", psMult:-ApproximatelyEqual(psMult, pusoConst), true);

Try[testnoerror]("test 14", PuSO:-BinaryMultiply(puso1, pusoConst), 'assign'='psMult');

pso := PowerSeries(10 + 15*(x + y) + 10*(x^2 + y^2) + 5*x*y^2);
Try[testnoerror]("test 15", PuSO(pso, [x=2,y=3]), 'assign'='checkMult');
Try("test 16", psMult:-ApproximatelyEqual(psMult, checkMult, 3), true);

Try[testnoerror]("test 17", PuSO:-BinaryMultiply(puso1, puso2), 'assign'='psMult');

pso := PowerSeries(4 + 6*(x^3 + y^5) + 4*(x^6 + y^(10)) + 2*x^3*y^(10)
				 + 4*x + 6*(x^4 + x*y^5) + 4*(x^7 + x*y^(10)) + 2*x^4*y^(10)
				 + 4*y^3 + 6*(x^3*y^3 + y^8) + 4*(x^6*y^3 + y^(13)) + 2*x^3*y^(13));
Try[testnoerror]("test 18", PuSO(pso, [x=6,y=15]), 'assign'='checkMult');
Try("test 19", psMult:-ApproximatelyEqual(psMult, checkMult, 4), true);

Try[testnoerror]("test 20", PuSO:-BinaryMultiply(puso3, puso4), 'assign'='psMult');

pso := PowerSeries((1/2 - (1/2)*x - (1/2)*y^1 + (1/2)*x^2 + x*y + (1/2)*y^2)
					*(1/2-(3/4)*x^2-(3/4)*y^6 +(5/8)*x^4 + (9/4)*x^2*y^6 +(5/8)*y^12));
Try[testnoerror]("test 21", PuSO(pso, [x=4,y=6]), 'assign'='checkMult');
Try("test 22", psMult:-ApproximatelyEqual(psMult, checkMult, 2), true);

#################################################################################
# Setup
pso := PowerSeries(2 + 3*(x + y) + 2*(x^2 + y^2) + x*y^2 + z^2);
expDivs := table([x=2, y=3]);
Try[testnoerror]("test 23", PuSO(pso, expDivs), 'assign'='puso5');

pso := PowerSeries(2 + 2*(x + y)+ x^2*w);
expDivs := table([x=6, y=5, w=2]);
Try[testnoerror]("test 24", PuSO(pso, expDivs), 'assign'='puso6');

pso := PowerSeries(d -> x^d/d!, analytic=exp(x));
expDivs := table([x=2]);
Try[testnoerror]("test 25", PuSO(pso, expDivs), 'assign'='puso7');

# Test Add
Try[testnoerror]("test 26", PuSO:-BinaryMultiply(puso5, puso6), 'assign'='psMult');
Try[verify,table]("test 27", psMult:-GetExpDivs(psMult), 
								table(':-sparse' = 1, [x=6, y=15, z=1, w=2]));

pso := PowerSeries(4 + 3*x^3 + 3*y^5 + 2*x^6 + 2*y^(10) + x^3*y^(10) + z ^2+ 2*x + 2*y^3 + x^2*w);
Try[testnoerror]("test 28", PuSO(pso, [x=6,y=15, w=2, z=1]), 'assign'='checkMult');
Try("test 29", psMult:-ApproximatelyEqual(psMult, checkMult), true);								

Try[testnoerror]("test 30", PuSO:-BinaryMultiply(puso6, puso7), 'assign'='psMult');

pso := PowerSeries((2+ 2*x+2*y+x^2*w)*(1+x^3+(1/2)*x^6+(1/6)*x^9));
Try[testnoerror]("test 31", PuSO(pso, [x=6,y=5, w=2]), 'assign'='checkMult');
Try("test 32", psMult:-ApproximatelyEqual(psMult, checkMult, 9), true);	

Try[testnoerror]("test 33", PuSO:-BinaryMultiply(puso3, puso7), 'assign'='psMult');

pso := PowerSeries((1/2 - (1/2)*x - (1/2)*y^1 + (1/2)*x^2 + x*y + (1/2)*y^2)*(1+ x^2+(1/2)*x^4));
Try[testnoerror]("test 34", PuSO(pso, [x=4,y=6]), 'assign'='checkMult');
Try("test 35", psMult:-ApproximatelyEqual(psMult, checkMult, 2), true);	

pso := PowerSeries(d -> x^d/d!, analytic=exp(x));

Try[testnoerror]("test 36", PuSO:-BinaryMultiply(puso7, pso), 'assign'='psMult');

Try[verify,table]("test 37", psMult:-GetExpDivs(psMult), 
								table(':-sparse' = 1, [x=2]));

pso := PowerSeries((1 + x^2 + (1/2)*x^4)*(1 + x + (1/2)*x^2));
Try[testnoerror]("test 38", PuSO(pso, [x=2]), 'assign'='checkMult');
Try("test 39", psMult:-ApproximatelyEqual(psMult, checkMult, 2), true);	

pso := PowerSeries(d -> x^d/d!, analytic=exp(x));

Try[testnoerror]("test 40", PuSO:-BinaryMultiply(pso, puso7), 'assign'='psMult');

Try[verify,table]("test 41", psMult:-GetExpDivs(psMult), 
								table(':-sparse' = 1, [x=2]));

pso := PowerSeries((1 + x^2 + (1/2)*x^4)*(1 + x + (1/2)*x^2));
Try[testnoerror]("test 42", PuSO(pso, [x=2]), 'assign'='checkMult');
Try("test 43", psMult:-ApproximatelyEqual(psMult, checkMult, 2), true);	

#################################################################################
# Setup

upops := UnivariatePolynomialOverPowerSeries([PowerSeries(1), PowerSeries(0), PowerSeries(x), 
										  PowerSeries(y), 1/PowerSeries(1 + x + y)], z);

Try[testnoerror]("test 45", PuSO:-BinaryMultiply(puso7, upops), 'assign'='psMult');	

pso := PowerSeries((1 + x^(1) + (1/2)*x^2 + (1/6)*x^3 +(1/24)*x^4)*(1+x^2*z^2 + y*z^3 + z^4));
Try[testnoerror]("test 46", PuSO(pso, [x=2, z=1]), 'assign'='checkmult');
Try("test 47", psMult:-ApproximatelyEqual(psMult, checkmult, 3), true);


Try[testnoerror]("test 48", PuSO:-BinaryMultiply(upops, puso7), 'assign'='psMult');	

Try("test 49", psMult:-ApproximatelyEqual(psMult, checkmult, 3), true);
#end test
