#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PuSO := MultivariatePowerSeries:-PuiseuxSeriesObject:
kernelopts(opaquemodules=true):

# Inverse
# Setup
ord := [x, y];
pso := PowerSeries(1/(1+u*v));
ordCV := [u, v];
e := [x=-5, y=3];
rays := [[1,0], [1,-1]];

Try[testnoerror]("test 1", PuSO(pso, ord, ordCV, rays, e), 'assign'='puso1');

pso := PowerSeries(1/(1+u));
mp := [u=x^(-1)*y^2];
e := [x=3, y=-4];

Try[testnoerror]("test 2", PuSO(pso, mp, e), 'assign'='puso2');

pso := PowerSeries(2 + 2*(u + v));
mp := [u=x^(-1)*y^2, v=y]; 
e := [x=3, y=2];

Try[testnoerror]("test 3", PuSO(pso, mp, e), 'assign'='puso3');

ord := [x1, x2, x3, x4, x5];
pso := PowerSeries(1/(1+u1*u2*u3*u4));
e := [x1=-5, x2=2];
rays := [[-2,1,2,-1,0], [3,2,1,2,3], [2,-2,-4,5,-1], [-3,4,0,-1,0]];
ordCV := [u1, u2, u3, u4];

Try[testnoerror]("test 21", PuSO(pso, ord, ordCV, rays), 'assign'='puso4');

# Test Inverse
Try[testnoerror]("test 22", PuSO:-Inverse(puso1), 'assign'='inv');
Try[testnoerror]("test 23", PuSO:-BinaryMultiply(puso1, inv), 'assign'='pusMult');
Try("test 24", PuSO:-Truncate(pusMult, 40), 1);

Try[testnoerror]("test 25", PuSO:-Inverse(puso2), 'assign'='inv');
Try[testnoerror]("test 26", PuSO:-BinaryMultiply(puso2, inv), 'assign'='pusMult');
Try("test 27", PuSO:-Truncate(pusMult, 40), 1);

Try[testnoerror]("test 28", PuSO:-Inverse(puso3), 'assign'='inv');
Try[testnoerror]("test 29", PuSO:-BinaryMultiply(puso3, inv), 'assign'='pusMult');
Try("test 30", PuSO:-Truncate(pusMult, 40), 1);

Try[testnoerror]("test 31", PuSO:-Inverse(puso4), 'assign'='inv');
Try[testnoerror]("test 32", PuSO:-BinaryMultiply(puso4, inv), 'assign'='pusMult');
Try("test 33", PuSO:-Truncate(pusMult, 40), 1);

# Setup PuSO without constant
ord := [x, y];
pso := PowerSeries(1/(1+u^2*v) -1);
ordCV := [u, v];
e := [x=-5, y=3];
rays := [[1,1], [2,-1]];

Try[testnoerror]("test 34", PuSO(pso, ord, ordCV, rays, e), 'assign'='puso5');

q := u/(u + v +u^2);
mp := [u=x, v=x*y^(-1)]; 
e := [x=3, y=2];

Try[testnoerror]("test 35", PuSO(q, mp, e), 'assign'='puso6');

q := u/(u + v + w +u^2*w);
mp := [u=x, v=x*y^(-1), w= x*y]; 
e := [x=3, y=-2];

Try[testnoerror]("test 36", PuSO(q, mp, e), 'assign'='puso7');

# Test Inverse
Try[testnoerror]("test 37", PuSO:-Inverse(puso5), 'assign'='inv');
Try[testnoerror]("test 38", PuSO:-BinaryMultiply(puso5, inv), 'assign'='pusMult');
Try("test 39", PuSO:-Truncate(pusMult, 40), 1);

#Test without the rational constructor
Try[testnoerror]("test 37.1", PuSO:-Inverse(puso5, useanalytic=false), 'assign'='inv');
Try[testnoerror]("test 38.1", PuSO:-BinaryMultiply(puso5, inv), 'assign'='pusMult');
Try("test 39", PuSO:-Truncate(pusMult, 40), 1);

Try[testnoerror]("test 40", PuSO:-Inverse(puso6), 'assign'='inv');
Try[testnoerror]("test 41", PuSO:-BinaryMultiply(puso6, inv), 'assign'='pusMult');
Try("test 42", PuSO:-Truncate(pusMult, 40), 1);

#Test without the rational constructor
Try[testnoerror]("test 40.1", PuSO:-Inverse(puso6, useanalytic=false), 'assign'='inv');
Try[testnoerror]("test 41.1", PuSO:-BinaryMultiply(puso6, inv), 'assign'='pusMult');
Try("test 42", PuSO:-Truncate(pusMult, 40), 1);

Try[testnoerror]("test 43", PuSO:-Inverse(puso7), 'assign'='inv');
Try[testnoerror]("test 44", PuSO:-BinaryMultiply(puso7, inv), 'assign'='pusMult');
Try("test 45", PuSO:-Truncate(pusMult, 40), 1);

#Test without the rational constructor
Try[testnoerror]("test 43.1", PuSO:-Inverse(puso7, useanalytic=false), 'assign'='inv');
Try[testnoerror]("test 44.1", PuSO:-BinaryMultiply(puso7, inv), 'assign'='pusMult');
Try("test 45", PuSO:-Truncate(pusMult, 40), 1);

# Setup PuSO undefined st its rays have positive weight
ord := [x, y];
bproc := proc(d) 
	if d mod 2 = 0 then
		return 0;
	else  
		return u^(d-1)*v;
	end if; 
end proc;  

pso := PowerSeries(bproc, variables={u,v});
ordCV := [u, v];
e := [x=-5, y=3];
rays := [[1,0], [9,-1]];

Try[testnoerror]("test 46", PuSO(pso, ord, ordCV, rays, e), 'assign'='puso8');

ord := [x, y];
bproc := proc(d) 
	if d = 0 then
		return 0;
	elif d = 1 then  
		return u + v;
	elif d = 3 then
		return u^2*v;
	else 
		return 0;
	end if; 
end proc;  

pso := PowerSeries(bproc, variables={u,v});
ordCV := [u, v];
e := [x=-2, y=3];
rays := [[1,0], [9,-1]];

Try[testnoerror]("test 47", PuSO(pso, ord, ordCV, rays, e), 'assign'='puso9');

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
rays := [[1,0], [9,-1]];

Try[testnoerror]("test 48", PuSO(pso, ord, ordCV, rays, e), 'assign'='puso10');

# Test inverse
Try[testnoerror]("test 49", PuSO:-Inverse(puso8), 'assign'='inv');
Try[testnoerror]("test 50", PuSO:-BinaryMultiply(puso8, inv), 'assign'='pusMult');
Try("test 51", PuSO:-Truncate(pusMult, 40), 1);

Try[testnoerror]("test 52", PuSO:-Inverse(puso9), 'assign'='inv');
Try[testnoerror]("test 53", PuSO:-BinaryMultiply(puso9, inv), 'assign'='pusMult');
Try("test 54", PuSO:-Truncate(pusMult, 40), 1);

Try[testnoerror]("test 55", PuSO:-Inverse(puso10), 'assign'='inv1');
Try[testnoerror]("test 56", PuSO:-BinaryMultiply(puso10, inv1), 'assign'='pusMult');
Try("test 57", PuSO:-Truncate(pusMult, 40), 1);

Try[testnoerror]("test 58", PuSO:-Inverse(inv1), 'assign'='inv');
Try[testnoerror]("test 59", PuSO:-BinaryMultiply(inv1, inv), 'assign'='pusMult');
Try("test 60", PuSO:-Truncate(pusMult, 40), 1);

# Setup PuSO undefined st its rays have a ZERO RAY
ord := [x, y];
bproc := proc(d) 
	if d mod 2 = 0 then
		return 0;
	else  
		return u^(d-1)*v;
	end if; 
end proc;  

pso := PowerSeries(bproc, variables={u,v});
ordCV := [u, v];
e := [x=-5, y=3];
rays := [[1,3], [1,-1]];

Try[testnoerror]("test 61", PuSO(pso, ord, ordCV, rays, e), 'assign'='puso8');

ord := [x, y];
bproc := proc(d) 
	if d = 0 then
		return 0;
	elif d = 1 then  
		return u + v;
	elif d = 3 then
		return u^2*v;
	else 
		return 0;
	end if; 
end proc;  

pso := PowerSeries(bproc, variables={u,v});
ordCV := [u, v];
e := [x=-2, y=3];
rays := [[2,1], [1,-1]];

Try[testnoerror]("test 62", PuSO(pso, ord, ordCV, rays, e), 'assign'='puso9');

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
rays := [[2,0], [1,-1]];

Try[testnoerror]("test 63", PuSO(pso, ord, ordCV, rays, e), 'assign'='puso10');

# Test inverse
Try[testnoerror]("test 64", PuSO:-Inverse(puso8), 'assign'='inv');
Try[testnoerror]("test 65", PuSO:-BinaryMultiply(puso8, inv), 'assign'='pusMult');
Try("test 66", PuSO:-Truncate(pusMult, 40), 1);

Try[testnoerror]("test 67", PuSO:-Inverse(puso9), 'assign'='inv');
Try[testnoerror]("test 68", PuSO:-BinaryMultiply(puso9, inv), 'assign'='pusMult');
Try("test 69", PuSO:-Truncate(pusMult, 40), 1);

Try[testnoerror]("test 70", PuSO:-Inverse(puso10), 'assign'='inv1');
Try[testnoerror]("test 71", PuSO:-BinaryMultiply(puso10, inv1), 'assign'='pusMult');
Try("test 71", PuSO:-Truncate(pusMult, 40), 1);

Try[testnoerror]("test 72", PuSO:-Inverse(inv1), 'assign'='inv');
Try[testnoerror]("test 73", PuSO:-BinaryMultiply(inv1, inv), 'assign'='pusMult');
Try("test 74", PuSO:-Truncate(pusMult, 40), 1);

# Test Errors in bounds
ord := [x, y, z];
pso := PowerSeries(v+u^30+w^20);
ordCV := [u, v,w];
rays := [[1,-1,0], [1,1,2], [1,2,-1]];

Try[testnoerror]("test 75", PuSO(pso, ord, ordCV, rays), 'assign'='puso11');

Try[testnoerror]("test 76", PuSO:-Inverse(puso11, useanalytic=false), 'assign'='inv');

Try[testnoerror]("test 77", PuSO:-BinaryMultiply(puso11, inv), 'assign'='pusMult');

Try[testerror]("test 78", PuSO:-Truncate(pusMult, 40));

Try[testnoerror]("test 79", PuSO:-Inverse(puso11, 10, 30, useanalytic=false), 'assign'='inv');
Try[testnoerror]("test 80", PuSO:-BinaryMultiply(puso11, inv), 'assign'='pusMult');
Try("test 81", PuSO:-Truncate(pusMult, 40), 1);

#########################
## Puiseux Series Test ##
#########################
#Inverse
# With constant.
ord := [x, y];
pso := PowerSeries(1/(1+u*v));
ordCV := [u, v];
e := [x=-5, y=3];
rays := [[1/10,0], [1/5,-1/5]];

Try[testnoerror]("test 82", PuSO(pso, ord, ordCV, rays, e), 'assign'='puso12');

Try[testnoerror]("test 83", PuSO:-Inverse(puso12), 'assign'='inv');
Try[testnoerror]("test 84", PuSO:-BinaryMultiply(puso12, inv), 'assign'='pusMult');
Try("test 85", PuSO:-Truncate(pusMult, 40), 1);

# Without constant - rational constructor.
q := u/(u + v + w +u^2*w);
mp := [u=x^(1/4), v=x^(3)*y^(2), w=x^(1/2)*y^(-1/4)]; 
e := [x=3, y=-2];

Try[testnoerror]("test 86", PuSO(q, mp, e), 'assign'='puso13');

Try[testnoerror]("test 87", PuSO:-Inverse(puso13), 'assign'='inv');
Try[testnoerror]("test 88", PuSO:-BinaryMultiply(puso13, inv), 'assign'='pusMult');
Try("test 89", PuSO:-Truncate(pusMult, 40), 1);

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

Try[testnoerror]("test 90", PuSO(pso, ord, ordCV, rays, e), 'assign'='puso14');

Try[testnoerror]("test 91", PuSO:-Inverse(puso14), 'assign'='inv1');
Try[testnoerror]("test 92", PuSO:-BinaryMultiply(puso14, inv1), 'assign'='pusMult');
Try("test 93", PuSO:-Truncate(pusMult, 40), 1);

Try[testnoerror]("test 94", PuSO:-Inverse(inv1), 'assign'='inv');
Try[testnoerror]("test 95", PuSO:-BinaryMultiply(inv1, inv), 'assign'='pusMult');
Try("test 96", PuSO:-Truncate(pusMult, 40), 1);

# Exp power series.
ord := [x, y];
pso := PowerSeries(d -> (u+v)^d/d!, analytic=exp(u+v));
ordCV := [u, v];
rays := [[1/4,0], [3/5,-2/5]];
e := [];

Try[testnoerror]("test 97", PuSO(pso, ord, ordCV, rays, e), 'assign'='pus_exp1');
Try[testnoerror]("test 98", PuSO:-Inverse(pus_exp1), 'assign'='inv');
Try[testnoerror]("test 99", PuSO:-BinaryMultiply(pus_exp1, inv), 'assign'='pusMult');
Try("test 100", PuSO:-Truncate(pusMult, 40), 1);

##########################
## Univariate inversion ##
##########################
ord := [x];
pso := PowerSeries(u+3*u^3-u^5);
ordCV := [u];
e := [x=1];
rays := [[1/2]];

Try[testnoerror]("test univariate 1", PuSO(pso, ord, ordCV, rays, e), 'assign'='pus1');
Try[testnoerror]("test univariate 2", PuSO:-Inverse(pus1), 'assign'='inv');
Try[testnoerror]("test univariate 3", PuSO:-BinaryMultiply(pus1, inv), 'assign'='pusMult');
Try("test univariate 4", PuSO:-Truncate(pusMult, 40), 1);

ord := [x];
bproc := proc(d) 
	if d = 0 then
		return 0;
	elif d = 30 then
		return u^(30);
	elif d = 40 then
		return u^(40);
	else 
		return 0;
	end if; 
end proc;  

pso := PowerSeries(bproc, variables={u});
ordCV := [u];
e := [x=-2];
rays := [[1/3]];

Try[testnoerror]("test univariate 5", PuSO(pso, ord, ordCV, rays, e), 'assign'='pus2');
Try[testerror]("test univariate 6", PuSO:-Inverse(pus2));
Try[testnoerror]("test univariate 7", PuSO:-Inverse(pus2, 40), 'assign'='inv');
Try[testnoerror]("test univariate 8", PuSO:-BinaryMultiply(pus2, inv), 'assign'='pusMult');
Try("test univariate 9", PuSO:-Truncate(pusMult, 40), 1);

##########################
##  Bivariate inversion ##
##########################
ord := [x, y];
pso := u^2*v^2-u^3*v^3 + u^20*v;
ordCV := [u, v];
e := [x=-5, y=3];
rays1 := [[1,1/2], [1/2,-1/3]];
rays2 := [[1,1/2], [1/3,-1/3]];

Try[testnoerror]("test bivariate 1", PuSO(pso, ord, ordCV, rays1, e), 'assign'='puso1');
Try[testnoerror]("test bivariate 3", PuSO:-Inverse(puso1), 'assign'='inv');
Try[testnoerror]("test bivariate 4", PuSO:-BinaryMultiply(puso1, inv), 'assign'='pusMult');
Try("test bivariate 5", PuSO:-Truncate(pusMult, 40), 1);

Try[testnoerror]("test bivariate 6", PuSO(pso, ord, ordCV, rays2, e), 'assign'='puso2');
Try[testnoerror]("test bivariate 8", PuSO:-Inverse(puso2), 'assign'='inv');
Try[testnoerror]("test bivariate 9", PuSO:-BinaryMultiply(puso2, inv), 'assign'='pusMult');
Try("test bivariate 10", PuSO:-Truncate(pusMult, 40), 1);

p:= 18*u^55*v^44+16*u^15*v^84+46*u^74*v^24-49*u^71*v^27-66*u^14*v^84-22*u^93*v^4+48*u^54*v^43-82*u^12*v^85-50*u^17*v^79-61*u^55*v^39+45*u^20*v^74-32*u^85*v^8+17*u^73*v^20+87*u^61*v^32-60*u^60*v^33+7*u^19*v^74-40*u^51*v^41-68*u^49*v^43-15*u^9*v^83-34*u^73*v^18+72*u^10*v^81+28*u^82*v^7+63*u^75*v^14+21*u^68*v^21-71*u^33*v^56-66*u^20*v^69+86*u^69*v^19+57*u^58*v^30-35*u^32*v^56+57*u^13*v^75-2*u^85*v^2-95*u^79*v^7+61*u^29*v^57-58*u^77*v^8-79*u^4*v^80+89*u^62*v^21+14*u^21*v^62+96*u^58*v^24-51*u^57*v^25+62*u^34*v^47-96*u^40*v^40-8*u^30*v^50-54*u^18*v^62-14*u^31*v^48+83*u^24*v^55-97*u^48*v^30+86*u^74*v^3+u^69*v^7-95*u^66*v^10+45*u^61*v^9+49*u^54*v^16-58*u^14*v^56+49*u^11*v^59+16*u^58*v^10-85*u^53*v^15-44*u^45*v^23-31*u^30*v^38-22*u^23*v^42-51*u^16*v^49-79*u^30*v^34-21*u^21*v^42-27*u^2*v^61+9*u^62-70*u^45*v^16+42*u^25*v^36-56*u^41*v^18-91*u^34*v^25-3*u^57*v-8*u^24*v^33+30*u^11*v^46-5*u^24*v^31+36*u^20*v^35-63*u^46*v^6+5*u^23*v^28-12*u^46*v^4+78*u^36*v^14+81*u^38*v^11+65*u^25*v^24-91*u^10*v^35+9*u^25*v^18+27*u^24*v^18+90*u^37*v^4+74*u^25*v^16+80*u^27*v^13-57*u^39+85*u^21*v^
18+5*u^8*v^30-36*u*v^37+37*u^22*v^14+58*u^34+29*u^8*v^26+97*u^13*v^18-67*u^7*v^24-51*u^11*v^18+88*u^6*v^23-26*u^12*v^9+5*u^7*v^11-29*u^13*v+68*u^5*v^6+95*u^4*v;

Try[testnoerror]("test bivariate 11", PuSO(p, ord, ordCV, rays2, e), 'assign'='puso3');
Try[testnoerror]("test bivariate 12", PuSO:-Inverse(puso3), 'assign'='inv');
Try[testnoerror]("test bivariate 13", PuSO:-BinaryMultiply(puso3, inv), 'assign'='pusMult');
Try("test bivariate 14", PuSO:-Truncate(pusMult, 100), 1);

##########################
## trivariate inversion ##
##########################
ord := [x, y, z];
p := -68*u^16*v^2*w^2+27*u^13*v^5*w^2-5*u^11*v^7*w^2+33*u^10*v*w^9-87*u^9*v^4*w^7+24*u^9*v^3*w^8-18*u^9*v*w^10-8*u^7*v^9*w^4+37*u^5*v^12*w^3-96*u^3*v^14*w^3+5*u^12*v^6*w+4*u^11*v^5*w^3-89*u^11*v*w^7-54*u^11*w^8-81*u^7*v^10*w^2+22*u^7*v^4*w^8-38*u^7*v^3*w^9-50*u^5*v^7*w^7-75*u^5*v^5*w^9+11*u^4*v^3*w^12+80*u^3*v*w^15+14*u^2*v^15*w^2-56*u^2*v^11*w^6+8*u^2*v^4*w^13+80*u^2*v^3*w^14-91*u*v^15*w^3-48*v^13*w^6+53*v^7*w^12+95*u^10*v^3*w^5-38*u^9*v*w^8+45*u^7*v^11-99*u^7*v^4*w^7-80*u^7*v*w^10+42*u^6*v^12+86*u^6*v^3*w^9+7*u^4*v^9*w^5-11*u^3*v^11*w^4+2*u*v^17-4*v^6*w^12+5*v^4*w^14-40*v^2*w^16-45*u^15*w^2-85*u^13*v*w^3+58*u^12*v^2*w^3-62*u^11*v^3*w^3-77*u^9*v^6*w^2-45*u^8*v*w^8-34*u^3*v^10*w^4+94*u^3*v^4*w^10+67*v^17+43*u^9*v^5*w^2+58*u^8*v^8-84*u^5*v^5*w^6+78*u^5*w^11-12*u^4*v^4*w^8-2*v^9*w^7+77*u^10*v*w^4-46*u^6*v^3*w^6+63*u^4*v*w^10+27*u^3*v^9*w^3-43*u^2*v^4*w^9+72*u^2*w^13-64*v^7*w^8+8*v^4*w^11-84*u^9*w^5-60*u^8*v^4*w^2-83*u^6*v^5*w^3+84*u*w^13-44*u^12*w+67*u^4*v*w^8-72*u^2*v^9*w^2+2*u*v^8*w^4+4*w^13+11*u^9*v*w^2-91*u^4*v^4*w^4-53*u^2*v^8*w^2+37*u^2*w^10+76*v^3*w^9+21*u^8*v^3+18*u^8*v^2*w+75*u^4*v^3*w^4+10*u^3*v^2*w^6+18*u*v^10+60*u*v^9*w-90*v^2*w^9+94*u^5*v^3*w^2+6*u^3*v^2*w^5+88*u^5*v^4+52*u^5*v^3+u^4*v^4+87*u^2*v^4*w^2+11*u^2*v^3*w^3-11*u^2*w^6-34*u^5*v^2-22*v^6-55*w^6-76*u^4*v+46*u^2*w-36*u*w^2;
ordCV := [u, v, w];
rays := [[1,-1,0], [1/3,1/2,2], [1/2,1,-1/2]];

Try[testnoerror]("test trivariate 1", PuSO(p, ord, ordCV, rays), 'assign'='puso1');
Try[testnoerror]("test trivariate 2", PuSO:-Inverse(puso1), 'assign'='inv');
Try[testnoerror]("test trivariate 3", PuSO:-BinaryMultiply(puso1, inv), 'assign'='pusMult');
Try("test trivariate 4", PuSO:-Truncate(pusMult, 100), 1);


########################
## Rational inversion ##
########################
mp := [u=x^(1/2), v=x^(1/3)*y^(-1/6), w=y^(1/2)*z^(-1/2)];
num:=52*u^2*w^2+10*u^3+40*w^3+54*v*w+12*w^2-98*u;
den := 11*u^4-47*u^3*v+39*u^3*w+50*u^2*v^2-91*u^2*v*w+41*u^2*w^2+25*u*v^3+99*u*v^2*w-8*u*v*w^2+36*u*w^3+20*v^4-91*v^3*w-35*v^2*w^2-29*v*w^3+31*w^4+99*u^3-26*u^2*v+91*u^2*w+48*u*v^2-9*u*v*w-37*u*w^2+17*v^3+29*v^2*w+66*v*w^2+55*w^3-76*u^2+86*u*v-90*u*w-17*v^2-72*v*w+6*w^2+47*u+76*v-96*w-15;


Try[testnoerror]("test ratfun 1", PuSO(num/den, mp), 'assign'='pus_rat');

Try[testnoerror]("test ratfun 2", PuSO:-Inverse(pus_rat), 'assign'='inv');
Try[testnoerror]("test ratfun 3", PuSO:-BinaryMultiply(pus_rat, inv), 'assign'='pusMult');
Try("test ratfun 4", PuSO:-Truncate(pusMult, 40), 1);

##
num := (u^2+u-6)*(52*u^2*w^2+10*u^3+40*w^3+54*v*w+12*w^2-98*u);
den2 := 2*(u^2+u-6)*den;

Try[testnoerror]("test ratfun 5", PuSO(num/den2, mp), 'assign'='pus_rat');

Try[testnoerror]("test ratfun 6", PuSO:-Inverse(pus_rat), 'assign'='inv');
Try[testnoerror]("test ratfun 7", PuSO:-BinaryMultiply(pus_rat, inv), 'assign'='pusMult');
Try("test ratfun 8", PuSO:-Truncate(pusMult, 40), 1);

##
num := 52*u^2*w^2+10*u^3+40*w^3+54*v*w+12*w^2-98*u;
den3 := 2*(u^2+u-6)*den;

Try[testnoerror]("test ratfun 9", PuSO(num/den3, mp), 'assign'='pus_rat');

Try[testnoerror]("test ratfun 10", PuSO:-Inverse(pus_rat), 'assign'='inv');
Try[testnoerror]("test ratfun 11", PuSO:-BinaryMultiply(pus_rat, inv), 'assign'='pusMult');
Try("test ratfun 12", PuSO:-Truncate(pusMult, 40), 1);

##
num := 52*u^2*w^2+10*u^3+40*w^3+54*v*w+12*w^2-98*u;
den3 := 2*sqrt(2)*den;

Try[testnoerror]("test ratfun 13", PuSO(num/den3, mp), 'assign'='pus_rat');

Try[testnoerror]("test ratfun 14", PuSO:-Inverse(pus_rat), 'assign'='inv');
Try[testnoerror]("test ratfun 15", PuSO:-BinaryMultiply(pus_rat, inv), 'assign'='pusMult');
Try("test ratfun 16", PuSO:-Truncate(pusMult, 40), 1);

##
num := u*Pi;
den3 := 2*(u+v+w);

Try[testnoerror]("test ratfun 17", PuSO(num/den3, mp), 'assign'='pus_rat');

Try[testnoerror]("test ratfun 18", PuSO:-Inverse(pus_rat), 'assign'='inv');
Try[testnoerror]("test ratfun 19", PuSO:-BinaryMultiply(pus_rat, inv), 'assign'='pusMult');
Try("test ratfun 20", PuSO:-Truncate(pusMult, 40), 1);

##
num := u;
den3 := 2*Pi*(u+v+w);

Try[testnoerror]("test ratfun 21", PuSO(num/den3, mp), 'assign'='pus_rat');

Try[testnoerror]("test ratfun 22", PuSO:-Inverse(pus_rat), 'assign'='inv');
Try[testnoerror]("test ratfun 23", PuSO:-BinaryMultiply(pus_rat, inv), 'assign'='pusMult');
Try("test ratfun 24", PuSO:-Truncate(pusMult, 40), 1);

##
num := u+v+w;
den3 := exp(2)*u^2 - exp(2);

Try[testnoerror]("test ratfun 25", PuSO(num/den3, mp), 'assign'='pus_rat');

Try[testnoerror]("test ratfun 26", PuSO:-Inverse(pus_rat), 'assign'='inv');
Try[testnoerror]("test ratfun 27", PuSO:-BinaryMultiply(pus_rat, inv), 'assign'='pusMult');
Try("test ratfun 28", PuSO:-Truncate(pusMult, 40), 1);
end test