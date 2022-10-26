#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
DeepCopy := MultivariatePowerSeries:-PuiseuxSeriesObject:-DeepCopy:
ConvertToUPoPS := MultivariatePowerSeries:-PuiseuxSeriesObject:-ConvertToUPoPS:
kernelopts(opaquemodules=true):

#Object definition.
pol := 2*(1/(1+y))*x^2 + 3*x*y + y;
var := 'x';

Try[testnoerror]("test obj def 1", UnivariatePolynomialOverPowerSeries(pol, var), 'assign'='upop1');
Try[testnoerror]("test obj def 2", Truncate(upop1, 5));

p2 := PowerSeries( 2*(1/(1+y)));
p1 := PowerSeries( 3*y);
p0 := PowerSeries( y);

Try[testnoerror]("test obj def 3", UnivariatePolynomialOverPowerSeries([p0,p1,p2], var), 'assign'='upop2');
Try[testnoerror]("test obj def 4", Truncate(upop2, 5));

s0 := PuiseuxSeries(p0);
A := Array([s0, p1, p2]);

Try[testnoerror]("test obj def 5", UnivariatePolynomialOverPowerSeries([s0,p1,p2], var), 'assign'='upop3');
Try[testnoerror]("test obj def 6", Truncate(upop3, 5));

Try[testnoerror]("test obj def 7", UnivariatePolynomialOverPowerSeries(A, var), 'assign'='upop4');

pol := 2*(1/(y+1))*x^2 + 3*x*y + y;
var := 'x';

Try[testnoerror]("test obj def 8", UnivariatePolynomialOverPowerSeries(pol, var), 'assign'='upop5');
Try[testnoerror]("test obj def 9", Truncate(upop5, 5));
Try("test obj def 10", upop5:-IsPuSOUPoP(upop5), false);

# Check error.
Try[testnoerror]("test obj def 11", WeierstrassPreparation(upop5));

Try[testnoerror]("test obj def 12", upop5:-ConvertToPuSOUPoP(upop5), 'assign'='upop5_puso');
Try[verify, normal]("test obj def 13", Truncate(upop5_puso, 5), Truncate(upop5,5));
Try("test obj def 14", upop5_puso:-IsPuSOUPoP(upop5_puso), true);

Try[testerror]("test obj def 15", WeierstrassPreparation(upop5_puso));

pol := 2*(1/y)*x^2 + 3*x*(1/y^2) + 1/(1+y);
var := 'x';

Try[testnoerror]("test obj def 16", UnivariatePolynomialOverPowerSeries(pol, var), 'assign'='upop6');
Try[testnoerror]("test obj def 17", Truncate(upop6, 5));
Try("test obj def 18", upop6:-IsPuSOUPoP(upop6), true);

pol := (3*x*y + y)/(x*y^2+y^3+x^4+7*x^4*y^7);
var := 'x';

Try[testerror]("test obj def 19", UnivariatePolynomialOverPowerSeries(pol, var), 'assign'='upop7');

pol := (3*x*y + y)/(y^2+y^3+7*y^7);
var := 'x';

Try[testnoerror]("test obj def 20", UnivariatePolynomialOverPowerSeries(pol, var), 'assign'='upop7');

Try("test obj def 21", upop5:-IsPuSOUPoP(upop7), true);

### ConvertToUPoPS.
puso := PuiseuxSeries((1/y)*x^2 + 3*x + y^2);
Try[testnoerror]("test convert 1", ConvertToUPoPS(puso, 'x'), 'assign'='puso_to_upop');
Try[verify,normal]("test convert  2", Truncate(puso_to_upop, 5), Truncate(puso,5));

## Setup.
Try[testnoerror]("test setup 1", UnivariatePolynomialOverPowerSeries(2*(1/(y+1))*x^2 + 3*x*y + y, 'x'), 'assign'='up1');
Try[testnoerror]("test setup 2", UnivariatePolynomialOverPowerSeries(13*x^10*z + 2*x*y + 2*x + 1, 'x'), 'assign'='up2');
Try[testnoerror]("test setup 3", UnivariatePolynomialOverPowerSeries(2*(1/y)*x^2 + 3*x*(1/y^2) + 1, 'x'), 'assign'='up3');
pso := PowerSeries( 2*(1/(1+y)));
puso1 := PuiseuxSeries(1/y);
puso2 := PuiseuxSeries((1/y)*x^2 + 3*x + y^2);
puso3 := PuiseuxSeries((1/y)*x^(2) + 3*x + y^2, [x=x^(1/2),y=y]);

### Addition.
Try[testnoerror]("test add 1", up1 + up3, 'assign'='up_sum');
Try[testnoerror]("test add 2", GetAnalyticExpression(up1) + GetAnalyticExpression(up3), 'assign'='ana');
Try[verify,normal]("test add 3", GetAnalyticExpression(up_sum), ana);

Try[testnoerror]("test add 4", up1 + up3 + pso, 'assign'='up_sum');
Try[testnoerror]("test add 5", GetAnalyticExpression(up1) + GetAnalyticExpression(up3)
		+ GetAnalyticExpression(pso), 'assign'='ana');
Try[verify,normal]("test add 6", GetAnalyticExpression(up_sum), ana);

Try[testnoerror]("test 7", up1 + up3 + pso + puso1, 'assign'='up_sum');
Try[testnoerror]("test 8", GetAnalyticExpression(up1) + GetAnalyticExpression(up3) 
        + GetAnalyticExpression(pso) + GetAnalyticExpression(puso1), 'assign'='ana');
Try[verify,normal]("test add 9", GetAnalyticExpression(up_sum), ana);

Try[testnoerror]("test add 10", up1 + pso + puso2, 'assign'='up_sum');
Try[testnoerror]("test add 11", GetAnalyticExpression(up1) + GetAnalyticExpression(pso)
		+ GetAnalyticExpression(puso2), 'assign'='ana');
Try[verify,normal]("test add 12", GetAnalyticExpression(up_sum), ana);

# Check error.
Try[testerror]("test add 13", up1+pso+puso3);

## Multiply.
up_prod := up1*up3;
Try[testnoerror]("test mult 1", GetAnalyticExpression(up1)*GetAnalyticExpression(up3), 'assign'='ana');
Try[verify,normal]("test Mult 2", GetAnalyticExpression(up_prod), ana);

up_prod := up1*up3*pso;
Try[testnoerror]("test mult 3", GetAnalyticExpression(up1)*GetAnalyticExpression(up3)*GetAnalyticExpression(pso), 'assign'='ana');
Try[verify,normal]("test Mult 4", GetAnalyticExpression(up_prod), ana);

up_prod := up1*up3*pso*puso1;
Try[testnoerror]("test mult 5", GetAnalyticExpression(up1)*GetAnalyticExpression(up3)*GetAnalyticExpression(pso)*GetAnalyticExpression(puso1), 'assign'='ana');
Try[verify,normal]("test Mult 6", GetAnalyticExpression(up_prod), ana);

up_prod := up1*pso*puso2;
Try[testnoerror]("test mult 7", GetAnalyticExpression(up1)*GetAnalyticExpression(pso)*GetAnalyticExpression(puso2), 'assign'='ana');
Try[verify,normal]("test Mult 8", GetAnalyticExpression(up_prod), ana);

# Check error.
Try[testerror]("test Mult 9", up1*pso*puso3);

## IsUnit
upop_zero := UnivariatePolynomialOverPowerSeries([PuiseuxSeries(0)], 'x');

bproc := proc(d) 
	if d  = 50 then
		return u^(50);
	else  
		return 0;
	end if; 
end proc; 
pso := PowerSeries(bproc, variables={u});

Try[testnoerror]("test IsUnit 1", PuiseuxSeries(pso), 'assign'='puso4');
upop_u := UnivariatePolynomialOverPowerSeries([puso4], 'x');

bproc2 := proc(d) 
	if d  = 20 then
		return 2*u^(20);
	else  
		return 0;
	end if; 
end proc; 
pso2 := PowerSeries(bproc2, variables={u});

Try[testnoerror]("test IsUnit 2", PuiseuxSeries(pso2), 'assign'='puso5');
upop_u2 := UnivariatePolynomialOverPowerSeries([puso5], 'x');

Try[testerror]("test IsUnit 3", upop6:-IsUnit(upop6));
Try("test IsUnit 4", upop5:-IsUnit(upop5), false);
Try[testerror]("test IsUnit 5", upop5_puso:-IsUnit(upop5_puso));
Try[testerror]("test IsUnit 6", upop_zero:-IsUnit(upop_zero));
Try[testerror]("test IsUnit 6", upop_u:-IsUnit(upop_u));
Try[testerror]("test IsUnit 7", upop_u2:-IsUnit(upop_u2));

## ApproximatelyZero
Try("test Appr zero 1", upop_u:-ApproximatelyZero(upop_u, 20), true);
Try("test Appr zero 2", upop_u:-ApproximatelyZero(upop_u, 60), false);
Try("test Appr zero 3", upop_zero:-ApproximatelyZero(upop_zero, 60), true);

## ApproximatelyEqual
## UPoPS and UPoPuS
Try("test Appr equal 1", upop_u:-ApproximatelyEqual(upop5, upop5_puso, 20), true);

## UPoPuS and UPoPus
Try("test Appr equal 2", upop_u:-ApproximatelyEqual(upop_zero, upop_u, 5), true);
Try("test Appr equal 3", upop_u:-ApproximatelyEqual(upop_zero, upop_u, 60), false);

##  GetCoefficient
Try("test get coeff 1", Truncate(upop5:-GetCoefficient(upop5, 1)), 3*y);

## Degree
Try("test Degree 1",upop5:-Degree(upop5), 2);

## UpdatePrecision
Try[testerror]("test UpdatePrecision 5", UpdatePrecision(upop5));

#end test