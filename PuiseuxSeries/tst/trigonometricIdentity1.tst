#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PuSO := MultivariatePowerSeries:-PuiseuxSeriesObject:
kernelopts(opaquemodules=true):

############
# Trigonometry power series.
ord := [x, y];
pso := PowerSeries(1);
ordCV := [];
rays := [];

Try[testnoerror]("test 1", PuSO(pso, ord, ordCV, rays), 'assign'='pusoOne');

ord := [x, y];
bproc1 := proc(d) 
	if d mod 2 = 0 then
		return (-1)^(d/2)*(u+v)^(d)/(d)!;
	else  
		return 0;
	end if; 
end proc;  

pso_cos := PowerSeries(bproc1, analytic=cos(u+v));
ordCV := [u, v];
rays := [[1/4,0], [3/5,-2/5]];
e := [];

Try[testnoerror]("test 2", PuSO(pso_cos, ord, ordCV, rays, e), 'assign'='pus_cos');
Try[testnoerror]("test 3", PuSO:-BinaryMultiply(pus_cos, pus_cos), 'assign'='pus_cos2');

ord := [x, y];
bproc2 := proc(d) 
	if d mod 2 = 0 then
		return 0;
	else  
		return (-1)^((d-1)/2)*(u+v)^(d)/(d)!;
	end if; 
end proc;  

pso_sin := PowerSeries(bproc2, analytic=sin(u+v));
ordCV := [u, v];
rays := [[1/4,0], [3/5,-2/5]];
e := [];

Try[testnoerror]("test 4", PuSO(pso_sin, ord, ordCV, rays, e), 'assign'='pus_sin');
Try[testnoerror]("test 5", PuSO:-BinaryMultiply(pus_sin, pus_sin), 'assign'='pus_sin2');

Try[testnoerror]("test 6", PuSO:-BinaryAdd(pus_cos2, pus_sin2), 'assign'='pusAdd');

Try("test 7", PuSO:-Truncate(pusAdd, 40), 1);

#end test