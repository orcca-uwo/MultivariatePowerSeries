#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PuSO := MultivariatePowerSeries:-PuiseuxSeriesObject:
kernelopts(opaquemodules=true):

############
ord := [x, y];
bproc1 := proc(d) 
	if d mod 2 = 0 then
		return (-1)^(d/2)*(u)^(d)/(d)!;
	else  
		return 0;
	end if; 
end proc;  

pso_cos := PowerSeries(bproc1, analytic=cos(u));
ordCV := [u];
rays := [[3/5,-2/5]];
e := [];

Try[testnoerror]("test 1", PuSO(pso_cos, ord, ordCV, rays, e), 'assign'='pus_cos');

ord := [x, y];
bproc2 := proc(d) 
	if d mod 2 = 0 then
		return 0;
	else  
		return (-1)^((d-1)/2)*(v)^(d)/(d)!;
	end if; 
end proc;  

pso_sin := PowerSeries(bproc2, analytic=sin(v));
ordCV := [v];
rays := [[3/5,-2/5]];
e := [];

Try[testnoerror]("test 2", PuSO(pso_sin, ord, ordCV, rays, e), 'assign'='pus_sin');

ord := [x, y];
bproc3 := proc(d) 
	if d mod 2 = 0 then
		return 0;
	else  
		return (-1)^((d-1)/2)*(u+v)^(d)/(d)!;
	end if; 
end proc;  

pso_sin := PowerSeries(bproc3, analytic=sin(u+v));
ordCV := [u,v];
rays := [[3/5,-2/5],[3/5,-2/5]];
e := [];

Try[testnoerror]("test 3", PuSO(pso_sin, ord, ordCV, rays, e), 'assign'='pus_sin2');

Try[testnoerror]("test 4", PuSO:-BinaryMultiply(pus_cos, pus_sin), 'assign'='pusMult');

Try[testnoerror]("test 5", PuSO:-BinaryAdd(pusMult, pusMult), 'assign'='pusAdd');

Try("test 6", pusMult:-ApproximatelyEqual(pusAdd, pus_sin2, 10, mode=absolute), true);

#end test