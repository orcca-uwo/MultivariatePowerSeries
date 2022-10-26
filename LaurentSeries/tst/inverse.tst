#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
LSO := MultivariatePowerSeries:-LaurentSeriesObject:
kernelopts(opaquemodules=true):

# Inverse
# Setup
ord := [x, y];
pso := PowerSeries(1/(1+u*v));
ordCV := [u, v];
e := [x=-5, y=3];
rays := [[1,0], [1,-1]];

Try[testnoerror]("test 1", LSO(pso, ord, ordCV, rays, e), 'assign'='lso1');

pso := PowerSeries(1/(1+u));
mp := [u=x^(-1)*y^2];
e := [x=3, y=-4];

Try[testnoerror]("test 2", LSO(pso, mp, e), 'assign'='lso2');

pso := PowerSeries(2 + 2*(u + v));
mp := [u=x^(-1)*y^2, v=y]; 
e := [x=3, y=2];

Try[testnoerror]("test 3", LSO(pso, mp, e), 'assign'='lso3');

ord := [x1, x2, x3, x4, x5];
pso := PowerSeries(1/(1+u1*u2*u3*u4));
e := [x1=-5, x2=2];
rays := [[-2,1,2,-1,0], [3,2,1,2,3], [2,-2,-4,5,-1], [-3,4,0,-1,0]];
ordCV := [u1, u2, u3, u4];

Try[testnoerror]("test 21", LSO(pso, ord, ordCV, rays), 'assign'='lso4');

# Test Inverse
Try[testnoerror]("test 22", LSO:-Inverse(lso1), 'assign'='inv');
Try[testnoerror]("test 23", LSO:-BinaryMultiply(lso1, inv), 'assign'='lsMult');
Try("test 24", LSO:-Truncate(lsMult, 40), 1);

Try[testnoerror]("test 25", LSO:-Inverse(lso2), 'assign'='inv');
Try[testnoerror]("test 26", LSO:-BinaryMultiply(lso2, inv), 'assign'='lsMult');
Try("test 27", LSO:-Truncate(lsMult, 40), 1);

Try[testnoerror]("test 28", LSO:-Inverse(lso3), 'assign'='inv');
Try[testnoerror]("test 29", LSO:-BinaryMultiply(lso3, inv), 'assign'='lsMult');
Try("test 30", LSO:-Truncate(lsMult, 40), 1);

Try[testnoerror]("test 31", LSO:-Inverse(lso4), 'assign'='inv');
Try[testnoerror]("test 32", LSO:-BinaryMultiply(lso4, inv), 'assign'='lsMult');
Try("test 33", LSO:-Truncate(lsMult, 40), 1);

# Setup LSO without constant
ord := [x, y];
pso := PowerSeries(1/(1+u^2*v) -1);
ordCV := [u, v];
e := [x=-5, y=3];
rays := [[1,1], [2,-1]];

Try[testnoerror]("test 34", LSO(pso, ord, ordCV, rays, e), 'assign'='lso5');

q := u/(u + v +u^2);
mp := [u=x, v=x*y^(-1)]; 
e := [x=3, y=2];

Try[testnoerror]("test 35", LSO(q, mp, e), 'assign'='lso6');

q := u/(u + v + w +u^2*w);
mp := [u=x, v=x*y^(-1), w= x*y]; 
e := [x=3, y=-2];

Try[testnoerror]("test 36", LSO(q, mp, e), 'assign'='lso7');

# Test Inverse
Try[testnoerror]("test 37", LSO:-Inverse(lso5), 'assign'='inv');
Try[testnoerror]("test 38", LSO:-BinaryMultiply(lso5, inv), 'assign'='lsMult');
Try("test 39", LSO:-Truncate(lsMult, 40), 1);

#Test without the rational constructor
Try[testnoerror]("test 37.1", LSO:-Inverse(lso5, useanalytic=false), 'assign'='inv');
Try[testnoerror]("test 38.1", LSO:-BinaryMultiply(lso5, inv), 'assign'='lsMult');
Try("test 39", LSO:-Truncate(lsMult, 40), 1);

Try[testnoerror]("test 40", LSO:-Inverse(lso6), 'assign'='inv');
Try[testnoerror]("test 41", LSO:-BinaryMultiply(lso6, inv), 'assign'='lsMult');
Try("test 42", LSO:-Truncate(lsMult, 40), 1);

#Test without the rational constructor
Try[testnoerror]("test 40.1", LSO:-Inverse(lso6, useanalytic=false), 'assign'='inv');
Try[testnoerror]("test 41.1", LSO:-BinaryMultiply(lso6, inv), 'assign'='lsMult');
Try("test 42", LSO:-Truncate(lsMult, 40), 1);

Try[testnoerror]("test 43", LSO:-Inverse(lso7), 'assign'='inv');
Try[testnoerror]("test 44", LSO:-BinaryMultiply(lso7, inv), 'assign'='lsMult');
Try("test 45", LSO:-Truncate(lsMult, 40), 1);

#Test without the rational constructor
Try[testnoerror]("test 43.1", LSO:-Inverse(lso7, useanalytic=false), 'assign'='inv');
Try[testnoerror]("test 44.1", LSO:-BinaryMultiply(lso7, inv), 'assign'='lsMult');
Try("test 45", LSO:-Truncate(lsMult, 40), 1);

# Setup LSO undefined st its rays have positive weight
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

Try[testnoerror]("test 46", LSO(pso, ord, ordCV, rays, e), 'assign'='lso8');

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

Try[testnoerror]("test 47", LSO(pso, ord, ordCV, rays, e), 'assign'='lso9');

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

Try[testnoerror]("test 48", LSO(pso, ord, ordCV, rays, e), 'assign'='lso10');

# Test inverse
Try[testnoerror]("test 49", LSO:-Inverse(lso8), 'assign'='inv');
Try[testnoerror]("test 50", LSO:-BinaryMultiply(lso8, inv), 'assign'='lsMult');
Try("test 51", LSO:-Truncate(lsMult, 40), 1);

Try[testnoerror]("test 52", LSO:-Inverse(lso9), 'assign'='inv');
Try[testnoerror]("test 53", LSO:-BinaryMultiply(lso9, inv), 'assign'='lsMult');
Try("test 54", LSO:-Truncate(lsMult, 40), 1);

Try[testnoerror]("test 55", LSO:-Inverse(lso10), 'assign'='inv1');
Try[testnoerror]("test 56", LSO:-BinaryMultiply(lso10, inv1), 'assign'='lsMult');
Try("test 57", LSO:-Truncate(lsMult, 40), 1);

Try[testnoerror]("test 58", LSO:-Inverse(inv1), 'assign'='inv');
Try[testnoerror]("test 59", LSO:-BinaryMultiply(inv1, inv), 'assign'='lsMult');
Try("test 60", LSO:-Truncate(lsMult, 40), 1);

# Setup LSO undefined st its rays have a ZERO RAY
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

Try[testnoerror]("test 61", LSO(pso, ord, ordCV, rays, e), 'assign'='lso8');

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

Try[testnoerror]("test 62", LSO(pso, ord, ordCV, rays, e), 'assign'='lso9');

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

Try[testnoerror]("test 63", LSO(pso, ord, ordCV, rays, e), 'assign'='lso10');

# Test inverse
Try[testnoerror]("test 64", LSO:-Inverse(lso8), 'assign'='inv');
Try[testnoerror]("test 65", LSO:-BinaryMultiply(lso8, inv), 'assign'='lsMult');
Try("test 66", LSO:-Truncate(lsMult, 40), 1);

Try[testnoerror]("test 67", LSO:-Inverse(lso9), 'assign'='inv');
Try[testnoerror]("test 68", LSO:-BinaryMultiply(lso9, inv), 'assign'='lsMult');
Try("test 69", LSO:-Truncate(lsMult, 40), 1);

Try[testnoerror]("test 70", LSO:-Inverse(lso10), 'assign'='inv1');
Try[testnoerror]("test 71", LSO:-BinaryMultiply(lso10, inv1), 'assign'='lsMult');
Try("test 71", LSO:-Truncate(lsMult, 40), 1);

Try[testnoerror]("test 72", LSO:-Inverse(inv1), 'assign'='inv');
Try[testnoerror]("test 73", LSO:-BinaryMultiply(inv1, inv), 'assign'='lsMult');
Try("test 74", LSO:-Truncate(lsMult, 40), 1);

# Test Errors in bounds
ord := [x, y, z];
pso := PowerSeries(v+u^30+w^20);
ordCV := [u, v,w];
rays := [[1,-1,0], [1,1,2], [1,2,-1]];

Try[testnoerror]("test 75", LSO(pso, ord, ordCV, rays), 'assign'='lso11');

Try[testnoerror]("test 76", LSO:-Inverse(lso11, useanalytic=false), 'assign'='inv');
Try[testnoerror]("test 77", LSO:-BinaryMultiply(lso11, inv), 'assign'='lsMult');
Try[testerror]("test 78", LSO:-Truncate(lsMult, 40));

Try[testnoerror]("test 79", LSO:-Inverse(lso11, 10, 30, useanalytic=false), 'assign'='inv');
Try[testnoerror]("test 80", LSO:-BinaryMultiply(lso11, inv), 'assign'='lsMult');
Try("test 81", LSO:-Truncate(lsMult, 40), 1);

#end test