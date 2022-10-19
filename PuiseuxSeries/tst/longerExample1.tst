#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PuSO := MultivariatePowerSeries:-PuiseuxSeriesObject:
AbsoluteDegreeBound := MultivariatePowerSeries:-PuiseuxSeriesObject:-AbsoluteDegreeBound:
MakeOrdCompatible := MultivariatePowerSeries:-PuiseuxSeriesObject:-MakeOrdCompatible:
ExtendPuiseuxSeriesObject := MultivariatePowerSeries:-PuiseuxSeriesObject:-ExtendPuiseuxSeriesObject:
FactorOutMonomial := MultivariatePowerSeries:-PowerSeriesObject:-FactorOutMonomial;
MakeRaysCompatible := MultivariatePowerSeries:-PuiseuxSeriesObject:-MakeRaysCompatible:
BinaryMakeCompatible := MultivariatePowerSeries:-PuiseuxSeriesObject:-BinaryMakeCompatible;
kernelopts(opaquemodules=true):

ord := [x, y];
pso := PowerSeries(1/(1+u*v));
ordCV := [u, v];
e := [x=-5, y=3];
rays := [[1/2,0], [1,-1/3]];

Try[testnoerror]("test 1", PuSO(pso, ord, ordCV, rays, e), 'assign'='puso1');

Try[testnoerror]("test 2", PuSO:-Inverse(puso1), 'assign'='puso2');

ord := [x, y];
pso := PowerSeries(2+u*v+u^2*v+u*v^2+v^3);
ordCV := [u, v];
e := [];
rays := [[1,0], [1,-1/2]];

Try[testnoerror]("test 3", PuSO(pso, ord, ordCV, rays, e), 'assign'='puso3');

Try[testnoerror]("test 4", PuSO:-BinaryMultiply(puso2, puso3), 'assign'='dummy1');

ord := [x, y];
pso := PowerSeries(1);
ordCV := [];
rays := [];

Try[testnoerror]("test 5", PuSO(pso, ord, ordCV, rays), 'assign'='pusoOne');

Try[testnoerror]("test 6", PuSO(-pso, ord, ordCV, rays), 'assign'='pusoNegOne');

Try[testnoerror]("test 7", PuSO:-BinaryMultiply(pusoNegOne, dummy1), 'assign'='dummy2');

Try[testnoerror]("test 8", PuSO:-BinaryAdd(pusoOne, dummy2), 'assign'='puso4');

Try[testnoerror]("test 9", PuSO:-Inverse(puso4), 'assign'='puso5');

Try[testnoerror]("test 10", PuSO:-BinaryMultiply(puso5, puso2), 'assign'='dummy3');

Try[testnoerror]("test 11", PuSO:-BinaryMultiply(pusoNegOne, puso3), 'assign'='dummy4');

Try[testnoerror]("test 12", PuSO:-BinaryAdd(puso1, dummy4), 'assign'='dummy5');

Try[testnoerror]("test 13", PuSO:-BinaryMultiply(dummy3, dummy5), 'assign'='puso_final');

Try("test 14", PuSO:-Truncate(puso_final, 50), 1);


#end test