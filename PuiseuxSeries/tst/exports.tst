#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
UPoPSObject := MultivariatePowerSeries:-UnivariatePolynomialOverPowerSeriesObject:
kernelopts(opaquemodules=true):

########################
## Object definitions ##
########################
q := u/(u^2*v +u^2*v^2);
mp := [u=x^(1/2), v=x^(1/3)*y^(-1/6)]; 
e := [x=3, y=2];

Try[testnoerror]("test obj def 1", PuiseuxSeries(q, mp, e), 'assign'='puso');

ana := x^(13/6)*y^(13/6)/(1+x^(1/3)*y^(-1/6));
Try[verify,normal]("test obj def 2", GetAnalyticExpression(puso), ana);

p := 1-u+u^2-u^3+u^4;
mp := [u= x^(1/3)*y^(-1/6)];
e := [x=13/6,y=13/6];
Try[testnoerror]("test obj def 3", PuiseuxSeries(p, mp, e), 'assign'='pusocheck');

Try("test obj def 4", ApproximatelyEqual(puso, pusocheck, 4, mode=absolute), true);

Try[testnoerror]("test obj def 5", PuiseuxSeries(u^2*v +u^2*v^2), 'assign'='puso1');

Try[testnoerror]("test obj def 6", PuiseuxSeries(PowerSeries(u^2*v +u^2*v^2)), 'assign'='puso2');

########################
## Truncate           ##
########################
ord := [x, y];
pso := PowerSeries(1/(1+u*v));
ordCV := [u, v];
e := [x=-5, y=3];
rays := [[1,0], [1,-1/2]];

Try[testnoerror]("test truncate 1", PuiseuxSeries(pso, ord, ordCV, rays, e), 'assign'='puso5');

poly := x^(-5)*y^3*(x^(12)*y^(-6/2)-x^(10)*y^(-5/2) + x^(8)*y^(-4/2) -x^(6)*y^(-3/2) +x^(4)*y^(-2/2) -x^2*y^(-1/2)+1 );

Try[verify,normal]("test truncate 2", Truncate(puso5, 9, mode=absolute), poly);

########################
## Equality           ##
########################
ord := [x, y];
bproc := proc(d::nonnegint, $) 
	if d mod 2 = 0 then
		return (-1)^(d/2)*(u*v)^(d/2)+(-1)^(d/2+1)*(u*v)^(d/2);
	else  
		return 0;
	end if; 
end proc;  

pso_tel := PowerSeries(bproc, variables={u,v});
ordCV := [u, v];
rays := [[1/4,0], [3/5,-2/5]];
e := [];

Try[testnoerror]("test equality 1", PuiseuxSeries(pso_tel, ord, ordCV, rays, e), 'assign'='pus_tel');

pus_tel;

Truncate(pus_tel, 9);

Try("test equality 2", ApproximatelyZero(pus_tel, 8), true);

########################
## Gets               ##
########################
q := u*(u*v+1)/(u*(1 +v));
mp := [u=x^(1/2), v=x*y^(-1/2)]; 
e := [x=3, y=2];

Try[testnoerror]("test gets 1", PuiseuxSeries(q, mp, e), 'assign'='puso1');

# Get internal power series.
Try[testnoerror]("test gets 2", GetPowerSeries(puso1), 'assign'='pso1');
Try("test gets 3", GetAnalyticExpression(pso1), (u*v+1)/(1+v));

# Get Puiseux series order.
Try("test gets 4", GetPuiseuxSeriesOrder(puso1), [x,y]);

# Get power series order.
Try("test gets 5", GetPowerSeriesOrder(puso1), [u,v]);

# Get the monomial that multiplies the PuSO.
Try("test gets 6", GetMonomial(puso1), x^3*y^2);

# Get the rays.
Try("test gets 7", GetRays(puso1), [[1/2,0], [1, -1/2]]);

# Get the change of variables applied to the internal ps.
Try("test gets 8", ChangeOfVariables(puso1), mp);

########################
## Addition           ##
########################
pso := PowerSeries(1/(1+u));
mp := [u=x^(-1/3)*y^2];
e := [x=3, y=-4];

Try[testnoerror]("test add 1", PuiseuxSeries(pso, mp, e), 'assign'='puso4');

pso := PowerSeries(2 + 2*(u + v));
mp := [u=x^(-1/2)*y, v=y];
e := [x=3, y=2];

Try[testnoerror]("test add 2", PuiseuxSeries(pso, mp, e), 'assign'='puso5');

Try[testnoerror]("test add 3", Add(puso4, puso5), 'assign'='pusAdd');

pso := 1 - u^2+u^4-u^6 + u^8 +2*v^6 +2*u^3*v^4+2*v^7;
mp := [u=x^(-1/6)*y, v=y];
e := [x=3, y=-4];
pusocheck := PuiseuxSeries(pso, mp, e);
Try("test add 4", ApproximatelyEqual(pusAdd, pusocheck, 7, mode=absolute), true);

## PuSO plus PSO
Try[testnoerror]("test add 5", Add(puso4, PowerSeries(1/(1+x+y))), 'assign'='pusAdd');

Try[testnoerror]("test add 6", puso5 + PowerSeries(1/(1+x*y+x^2)), 'assign'='pusAdd');

## PuSO plus PSO plus algebraic
Try[testnoerror]("test add 7", puso5 + PowerSeries(1/(1+x*y+x^2)) + x
							   + (x*y+5*x)/(x+y), 'assign'='pusAdd');

Try[testnoerror]("test add 8", puso5 + PowerSeries(1/(1+x*y+x^2)) 
							   + (x*y+x^2)/(1+x*y+x^2), 'assign'='pusAdd');

Try[testnoerror]("test add 9", PowerSeries(1/(1+x*y+x^2)) + puso5, 'assign'='pusAdd');

## UPoPS plus PSO and PuSO
upop1 := UnivariatePolynomialOverPowerSeries(13*x^10 + 2*x*y + 2*x + 1, 'x');
pso := PowerSeries(1/(1+y*z+y^2));

Try[testnoerror]("test add 10", upop1 +pso);

Try[testnoerror]("test add 11", pso + upop1);

pso := 1 - u^2+u^4-u^6 + u^8 +2*v^6 +2*u^3*v^4+2*v^7;
puso7 := PuiseuxSeries(pso, [u=x*y, v=y], [x=3, y=-4]);

Try[testnoerror]("test add 12", upop1+ puso7);

Try[testnoerror]("test add 13", puso7+ upop1);


########################
## Multiplication     ##
########################
ord := [x, y];
pso := PowerSeries(1/(1+u*v));
ordCV := [u, v];
e := [x=-5, y=3];
rays := [[1/2,0], [1,-1/3]];

Try[testnoerror]("test multiply 1", PuiseuxSeries(pso, ord, ordCV, rays, e), 'assign'='puso4');

ord := [x, y];
pso := PowerSeries(2+u*v+u^2*v+u*v^2+v^3);
ordCV := [u, v];
e := [];
rays := [[1,0], [1,-1/2]];

Try[testnoerror]("test multiply 2", PuiseuxSeries(pso, ord, ordCV, rays, e), 'assign'='puso5');

ord := [x, y];
pso := PowerSeries(1+u*v+v);
ordCV := [u, v];
e := [];
rays := [[1/2,0], [1,-1/3]];

Try[testnoerror]("test multiply 3", PuiseuxSeries(pso, ord, ordCV, rays, e), 'assign'='puso6');

Try[testnoerror]("test multiply 4", Multiply(puso5, puso6), 'assign'='pusMult');

p := normal(eval(2+u*v+u^2*v+u*v^2+v^3, [u=x,v=x*y^(-1/2)]));
q :=normal(eval(1+u*v+v, [u=x^(1/2), v=x*y^(-1/3)]));

Try[verify,normal]("test multiply 5", normal(p*q- Truncate(pusMult, 50)), 0);

Try[testnoerror]("test multiply 6", Multiply(puso4, puso6), 'assign'='pusMult');

poly := (1-(u*v)+ (u*v)^2 -(u*v)^3+ (u*v)^4 -(u*v)^5)*(1+u*v+v);
mp := [u=x^(1/2), v=x*y^(-1/3)];
e := [x=-5, y=3];
pusocheck := PuiseuxSeries(poly, mp, e);
Try("test multiply 7", ApproximatelyEqual(pusMult, pusocheck, 4, mode=absolute), true);

Try[testnoerror]("test multiply 8", Multiply(puso4, puso4), 'assign'='pusMult');

poly := (1-(u*v)+ (u*v)^2 -(u*v)^3+ (u*v)^4 -(u*v)^5)^2;
mp := [u=x^(1/2), v=x*y^(-1/3)];
e := [x=-10, y=6];
pusocheck := PuiseuxSeries(poly, mp, e);
Try("test multiply 9", ApproximatelyEqual(pusMult, pusocheck, 4, mode=absolute), true);

# Exp power series.
ord := [x, y];
pso := PowerSeries(d -> (u+v)^d/d!, analytic=exp(u+v));
ordCV := [u, v];
rays := [[1/4,0], [3/5,-2/5]];
e := [];

Try[testnoerror]("test multiply 10", PuiseuxSeries(pso, ord, ordCV, rays, e), 'assign'='pus_exp1');

ord := [x, y];
pso := PowerSeries(d -> (v)^d/d!, analytic=exp(v));
ordCV := [v];
rays := [[3/5,-2/5]];
e := [];

Try[testnoerror]("test multiply 11", PuiseuxSeries(pso, ord, ordCV, rays, e), 'assign'='pus_exp2');

ord := [x, y];
pso := PowerSeries(d -> (u)^d/d!, analytic=exp(u));
ordCV := [u];
rays := [[1/4,0]];
e := [];

Try[testnoerror]("test multiply 12", PuiseuxSeries(pso, ord, ordCV, rays, e), 'assign'='pus_exp3');

Try[testnoerror]("test multiply 13", Multiply(pus_exp2, pus_exp3), 'assign'='pusMult');

Try("test multiply 14", ApproximatelyEqual(pusMult, pus_exp1, 5, mode=absolute), true);

## PSO times PuSO
Try[testnoerror]("test multiply 13", Multiply(pus_exp2, PuiseuxSeries(PowerSeries(1/(1+x+y)))), 'assign'='pusMult');

pso := PowerSeries(1/(1+x^2*y^3));
pso1 := PowerSeries(1+x*y);
cv := [x=x^(1/2), y=y^(1/3)];

puso1 := PuiseuxSeries(pso, cv);

Try[testnoerror]("test multiply 14", puso1*pso1, 'assign'='pusMult');
Try[verify,normal]("test multiply 15", Truncate(pusMult, 20), 1);

## UPoPS times PSO times PuSO
upop1 := UnivariatePolynomialOverPowerSeries(13*x^10 + 2*x*y + 2*x + 1, 'x');
pso := PowerSeries(1/(1+y*z+y^2));

Try[testnoerror]("test multiply 16", upop1*pso);

Try[testnoerror]("test multiply 17", pso*upop1);

Try[testnoerror]("test multiply 18", upop1*puso7);

Try[testnoerror]("test multiply 19", puso7*upop1);

upop1 := UnivariatePolynomialOverPowerSeries(13*x^10 + 2*x*y + 2*x + 1, 'x');

ord := [x, y];
pso := PowerSeries(2+u*v+u^2*v+u*v^2+v^3);
ordCV := [u, v];
e := [x=-5,y=3/2];
rays := [[1/3,0], [1,-1/2]];

puso1 := PuiseuxSeries(pso, ord, ordCV, rays, e);
puso2 := PuiseuxSeries(13*x^10 + 2*x*y + 2*x + 1);
pso2 := PowerSeries(13*x^10 + 2*x*y + 2*x + 1);

# upop1*puso1 not covered yet.
Try[testerror]("test multiply 20", Multiply(upop1,puso1), 'assign'='pusMul1');
Try[testnoerror]("test multiply 21", Multiply(puso2,puso1), 'assign'='pusMul2');
Try[testnoerror]("test multiply 22", Multiply(pso2,puso1), 'assign'='pusMul3');

#Try[verify,normal]("test extra add 23", GetAnalyticExpression(pusMul1), 
#										GetAnalyticExpression(pusMul2));
Try[verify,normal]("test extra add 24", GetAnalyticExpression(pusMul2), 
										GetAnalyticExpression(pusMul3));
#Try[verify,normal]("test extra add 25", GetAnalyticExpression(pusMul3), 
#										GetAnalyticExpression(pusMul1));	
########################
## Inverse            ##
########################
# With constant.
ord := [x, y];
pso := PowerSeries(1/(1+u*v));
ordCV := [u, v];
e := [x=-5, y=3];
rays := [[1/10,0], [1/5,-1/5]];

Try[testnoerror]("test inverse 1", PuiseuxSeries(pso, ord, ordCV, rays, e), 'assign'='puso12');

Try[testnoerror]("test inverse 2", Inverse(puso12), 'assign'='inv');
Try[testnoerror]("test inverse 3", Multiply(puso12, inv), 'assign'='pusMult');
Try("test inverse 4", Truncate(pusMult, 40), 1);

# Without constant - rational constructor.
q := u/(u + v + w +u^2*w);
mp := [u=x^(1/4), v=x^(3)*y^(2), w=x^(1/2)*y^(-1/4)]; 
e := [x=3, y=-2];

Try[testnoerror]("test inverse 5", PuiseuxSeries(q, mp, e), 'assign'='puso13');

Try[testnoerror]("test inverse 6", Inverse(puso13), 'assign'='inv');
Try[testnoerror]("test inverse 7", Multiply(puso13, inv), 'assign'='pusMult');
Try("test inverse 8", Truncate(pusMult, 40), 1);

#Test without the rational constructor
ord := [x, y];
bproc := proc(d) 
	if d = 0 then
		return 0;
	elif d = 3 then
		return u^2*v + u*v^2;
	elif d = 4 then
		return u^4;
	else 
		return 0;
	end if; 
end proc;  

pso := PowerSeries(bproc, variables={u,v});
ordCV := [u, v];
e := [x=-2, y=3];
rays := [[1/15,0], [2,-1/10]];

Try[testnoerror]("test inverse 9", PuiseuxSeries(pso, ord, ordCV, rays, e), 'assign'='puso14');

Try[testnoerror]("test inverse 10", Inverse(puso14), 'assign'='inv1');
Try[testnoerror]("test inverse 11", Multiply(puso14, inv1), 'assign'='pusMult');
Try("test inverse 12", Truncate(pusMult, 40), 1);

Try[testnoerror]("test inverse 13", Inverse(inv1), 'assign'='inv');
Try[testnoerror]("test inverse 14", Multiply(inv1, inv), 'assign'='pusMult');
Try("test inverse 15", Truncate(pusMult, 40), 1);

# Exp power series.
ord := [x, y];
pso := PowerSeries(d -> (u+v)^d/d!, analytic=exp(u+v));
ordCV := [u, v];
rays := [[1/4,0], [3/5,-2/5]];
e := [];

Try[testnoerror]("test inverse 16", PuiseuxSeries(pso, ord, ordCV, rays, e), 'assign'='pus_exp1');
Try[testnoerror]("test inverse 17", Inverse(pus_exp1), 'assign'='inv');
Try[testnoerror]("test inverse 18", Multiply(pus_exp1, inv), 'assign'='pusMult');
Try("test inverse 19", Truncate(pusMult, 40), 1);

Try[testerror]("test inverse 20", Inverse(pso, useanalytic=false), 'assign'='inv');
Try[testerror]("test inverse 21", Inverse(pso, 12), 'assign'='inv');

###################
##	UPOPS to PSO ##
###################
upop1 := UnivariatePolynomialOverPowerSeries(13*x^10 + 2*x*y + 2*x + 1, 'x');
upop1_as_pso := upop1:-ConvertToPowerSeries(upop1);

Try("test UPOPS to PSO 1", Variables(upop1_as_pso), {x,y});
#end test