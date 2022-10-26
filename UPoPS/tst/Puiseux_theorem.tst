#test 500
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
FactorVariable := MultivariatePowerSeries:-UnivariatePolynomialOverPowerSeriesObject:-FactorVariable;
kernelopts(opaquemodules=true):

##
test_pf := proc(f, d, mode:=automatic, bnd::nonnegint := FAIL, $) 
	local up, res, differ; 
	if bnd<> FAIL then
		up := PuiseuxTheorem(f, bnd, ':-returnleadingcoefficient'=mode); 
	else 
		up := PuiseuxTheorem(f, ':-returnleadingcoefficient'=mode); 
	end if;

	res := Multiply(seq(up));
	differ := res + ((-1)*f);
	res := Truncate(differ, d);
	#res := simplify(res);
	res := evala(res); 

	if res = 0 then 
		return true; 
	else 
		return res; 
	end if; 
end proc:

## TschirnhausenTransformation
f:= 1+ (3/y^2)*x + x^2;
up1 := UnivariatePolynomialOverPowerSeries(f, 'x');
a := GetCoefficient(up1, 1);
Try[testnoerror]("test Tschirnhausen 1", up1:-TschirnhausenTransformation(up1, a), 'assign'='transf');
Try[testnoerror]("test Tschirnhausen 2", transf:-GetCoefficient(transf, 1), 'assign'='coef');
Try("test Tschirnhausen 3", GetAnalyticExpression(coef), 0);

f_transf := eval(f, [x=x-(3/(2*y^2))]);
Try[verify,normal]("test Tschirnhausen 4", GetAnalyticExpression(transf), f_transf);

Try[testnoerror]("test Tschirnhausen 5", transf:-TschirnhausenTransformation(transf, detransform=true, a), 'assign'='detransf');
Try[verify,normal]("test Tschirnhausen 6", GetAnalyticExpression(detransf), 						  GetAnalyticExpression(up1));

up2 := UnivariatePolynomialOverPowerSeries(x^10 + 2*x*y + 2*x + 1, 'x');
a := GetCoefficient(up2, 9);
Try[testnoerror]("test Tschirnhausen 7", up2:-TschirnhausenTransformation(up2, a), 'assign'='transf');
Try[testnoerror]("test Tschirnhausen 8", transf:-GetCoefficient(transf, 9), 'assign'='coef');
Try("test Tschirnhausen 9", GetAnalyticExpression(coef), 0);

f:= x^10 + (1/y^3)*x^9 -2*x*y + 2*x + y^2;
up3 := UnivariatePolynomialOverPowerSeries(f, 'x');
a := GetCoefficient(up3, 9);
Try[testnoerror]("test Tschirnhausen 10", up3:-TschirnhausenTransformation(up3, a), 'assign'='transf');
Try[testnoerror]("test Tschirnhausen 11", transf:-GetCoefficient(transf, 9), 'assign'='coef');
Try("test Tschirnhausen 12", GetAnalyticExpression(coef), 0);

f_transf := eval(f, [x=x-(1/(10*y^3))]);
Try[verify,normal]("test Tschirnhausen 13", GetAnalyticExpression(transf), f_transf);

Try[testnoerror]("test Tschirnhausen 14", transf:-TschirnhausenTransformation(transf, detransform=true, a), 'assign'='detransf');
Try[verify,normal]("test Tschirnhausen 15", GetAnalyticExpression(detransf),						  GetAnalyticExpression(up3));

f := (y-1)/(y^3+y^2+y)*x^10 + 1/(y-1)*x^9 -2*x*y + 2*x + 1;
up4 := UnivariatePolynomialOverPowerSeries(f, 'x');
a := GetCoefficient(up4, 9);
Try("test Tschirnhausen 16", up4:-IsMonic(up4), false);
Try[testnoerror]("test Tschirnhausen 17", up4:-TschirnhausenTransformation(up4, a));

up4 := up4:-MakeMonic(up4);
f:= GetAnalyticExpression(up4);
a9 := up4:-GetCoefficient(up4, 9);
Try[testnoerror]("test Tschirnhausen 18", up4:-TschirnhausenTransformation(up4, a9), 'assign'='transf');
a9 := GetAnalyticExpression(a9);
Try[testnoerror]("test Tschirnhausen 19", transf:-GetCoefficient(transf, 9), 'assign'='coef');
Try("test Tschirnhausen 20", GetAnalyticExpression(coef), 0);

f_transf := eval(f, [x=x-(1/10)*a9]);
Try[verify,normal]("test Tschirnhausen 21", GetAnalyticExpression(transf), f_transf);

# Ord of a PuSO.
puso1 := PuiseuxSeries(1/(1-u), [u=x^(1/3)], [x=-6]);
Try("test order 1", puso1:-GetOrder(puso1), -6);

puso2 := PuiseuxSeries(1/(1-u^7)-1, [u=x^(1/3)], [x=-6]);
Try("test order 2", puso2:-GetOrder(puso2), -11/3);

# get error.
puso3 := PuiseuxSeries(u^(60), [u=x^(1/3)], [x=-6]);
Try[testerror]("test order 3", puso3:-GetOrder(puso3));
Try("test order 4", puso3:-GetOrder(puso3, 60), 14);

## Puiseux tehorem
print("Puiseux theorem");
Try[testnoerror]("test Puiseux fac 1", up1:-PuiseuxTheorem(up1));

p2 := PowerSeries( 1/(1+y));
p1 := PowerSeries( 3*y);
s0 := PuiseuxSeries(PowerSeries(y), [y=y^(1/3)], [y=-5]);
A := Array([s0, p1, p2]);

Try[testnoerror]("test Puiseux fac 2", UnivariatePolynomialOverPowerSeries(A, 'x'), 'assign'='up5');

Try[testnoerror]("test Puiseux fac 3", up3:-PuiseuxTheorem(up3, 10));
Try[testnoerror]("test Puiseux fac 4", up5:-PuiseuxTheorem(up5, 10));

## Correctness of the factorization
p2 := PowerSeries( 1/(1+y));
p1 := PowerSeries( 3*y);
p0 := PowerSeries( y);

s0 := PuiseuxSeries(p0, [y=y^(1/3)], [y=-5]);
A := Array([s0, p1, PowerSeries(1)]);

Try[testnoerror]("test Puiseux theo 1", UnivariatePolynomialOverPowerSeries(A, 'x'), 'assign'='up5');

Try[testnoerror]("test Puiseux theo 2", UnivariatePolynomialOverPowerSeries(x^4-y^2, 'x'), 'assign'='up6');

Try[testnoerror]("test Puiseux theo 3", UnivariatePolynomialOverPowerSeries(x^3+3*x^2-y^2*x-3*y^2, 'x'), 'assign'='up7');

Try[testnoerror]("test Puiseux theo 4", UnivariatePolynomialOverPowerSeries([s0,p1,p2], 'x'), 'assign'='up8');

Try[testnoerror]("test Puiseux theo 5", UnivariatePolynomialOverPowerSeries(x^4+(y/(1+y^2))*x^2+y*x+1/y^2, 'x'), 'assign'='up9');

a0 := PuiseuxSeries((1+y^2)/y, [y=y^(1/3)]);
a1 := PuiseuxSeries(1 + y^2, [y=y^(1/2)], [y=-2]);
a2 := PuiseuxSeries(1/(y^2+1), [y=y^(1/6)], [y=-1]);
a3 := PuiseuxSeries(1/(y^2+y+1), [y=y^(1/2)], [y=-2]);

Try[testnoerror]("test Puiseux theo 6", UnivariatePolynomialOverPowerSeries([a0,a1,a2,a3], 'x'), 'assign'='up10');

Try[testnoerror]("test Puiseux theo 7", UnivariatePolynomialOverPowerSeries(x^3-y*x+y, 'x'), 'assign'='up11');

Try[testnoerror]("test Puiseux theo 8", UnivariatePolynomialOverPowerSeries(x^3 - x^2 * y^2-x*y^3 + y^4, 'x'), 'assign'='up12');

Try[testnoerror]("test Puiseux theo 9", UnivariatePolynomialOverPowerSeries(y^2*(x+2)-x^2*(2-x), 'y'), 'assign'='up13');

p:= (y^2 - x^2)*(x-1)*(2*x-3) - 4*(x^2 + y^2 - 2*x)^2;
Try[testnoerror]("test Puiseux theo 10", UnivariatePolynomialOverPowerSeries(p, 'y'), 'assign'='up14');

p:= (x^2 + y^2 - 1)^3 + 27*x^2*y^2;
Try[testnoerror]("test Puiseux theo 11", UnivariatePolynomialOverPowerSeries(p, 'y'), 'assign'='up15');

Try("test Puiseux theo 12", test_pf(up5, 20), true);

Try("test Puiseux theo 13", test_pf(up6, 20), true);

Try("test Puiseux theo 14", test_pf(up7, 20), true);

Try("test Puiseux theo 15", test_pf(up8, 20, true), true);

Try("test Puiseux theo 15.1", test_pf(up9, 10), true);

Try("test Puiseux theo 16", test_pf(up10, 20, true), true);

Try("test Puiseux theo 17", test_pf(up11, 20), true);

Try("test Puiseux theo 18", test_pf(up12, 20), true);

Try("test Puiseux theo 19", test_pf(up13, 20, true), true);

Try("test Puiseux theo 20", test_pf(up14, 20, true), true);

Try("test Puiseux theo 21", test_pf(up15, 20), true);

# Factor variable test.
a0 := PuiseuxSeries((1+y^2)/y, [y=y^(1/3)]);
a1 := PuiseuxSeries(1 + y^2, [y=y^(1/2)], [y=-2]);
a2 := PuiseuxSeries(1/(y^2+1), [y=y^(1/6)], [y=-1]);
a3 := PuiseuxSeries(1/(y^2+y+1), [y=y^(1/2)], [y=-2]);

Try[testnoerror]("test Factor Variable 1", UnivariatePolynomialOverPowerSeries(x^3-y*x^2+y*x, 'x'), 'assign'='upo1');
Try[testnoerror]("test Factor Variable 2", UnivariatePolynomialOverPowerSeries(x^4-y*x^3+y*x^2, 'x'), 'assign'='upo2');

bproc := proc(d) 
	if d  = 20 then
		return u^(20);
	else  
		return 0;
	end if; 
end proc; 
pso := PowerSeries(bproc, variables={u});
a0 := PuiseuxSeries(pso, [u=y]);
Try[testnoerror]("test Factor Variable 3", UnivariatePolynomialOverPowerSeries([a0,a1,a2,a3], 'x'), 'assign'='upo3');
Try[testnoerror]("test Factor Variable 4", UnivariatePolynomialOverPowerSeries([PuiseuxSeries(0),a1,a2,a3], 'x'), 'assign'='upo4');
Try[testnoerror]("test Factor Variable 5", UnivariatePolynomialOverPowerSeries([PuiseuxSeries(0), PuiseuxSeries(0),a0,a2,a3], 'x'), 'assign'='upo5');

Try[testnoerror]("test Factor Variable 6", FactorVariable(upo1), 'assign'='facs');
Try("test Factor Variable 7", facs[1], 1);
Try("test Factor Variable 8", Degree(facs[2]), 2);

Try[testnoerror]("test Factor Variable 9", FactorVariable(upo2), 'assign'='facs');
Try("test Factor Variable 10", facs[1], 2);
Try("test Factor Variable 11", Degree(facs[2]), 2);

Try[testnoerror]("test Factor Variable 12", FactorVariable(upo4), 'assign'='facs');
Try("test Factor Variable 13", facs[1], 1);
Try("test Factor Variable 14", Degree(facs[2]), 2);

Try[testerror]("test Factor Variable 15", FactorVariable(upo5));

Try[testnoerror]("test Factor Variable 26", FactorVariable(upo5, 21), 'assign'='facs');
Try("test Factor Variable 17", facs[1], 2);
Try("test Factor Variable 18", Degree(facs[2]), 2);

Try("test Factor Variable 19", test_pf(upo1, 20), true);
Try("test Factor Variable 20", test_pf(upo2, 20), true);
Try("test Factor Variable 21", test_pf(upo3, 20, true, 21), true);
Try("test Factor Variable 22", test_pf(upo4, 20, true), true);
Try("test Factor Variable 23", test_pf(upo5, 30, true, 21), true);

#end test