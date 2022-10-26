#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

# Setup.
puso1 := PuiseuxSeries(5);
puso2 := PuiseuxSeries(1/(1-u), [u=x^3], [x=-6]);
puso3 := PuiseuxSeries(u^2+3*u-1, [u=x^2]);
puso4 := PuiseuxSeries(u^2*v+u*v+u+v, [u=x^(2),v=x*y], [x=2]);
puso5 := PuiseuxSeries((u^4+3*u^2-1)/(1-u^6), [u=x^(1/2)]);
puso6 := PuiseuxSeries(1/(1-u), [u=x^3]);

bproc := proc(d) 
	if d mod 2 = 0 then
		return 0;
	else  
		return u^(d-1)*v;
	end if; 
end proc;  

pso := PowerSeries(bproc, variables={u,v});
puso7 := PuiseuxSeries(pso, [x,y], [u,v], [[0,1], [1,-1]], [x=-5, y=3]);
puso8 := PuiseuxSeries(pso, [x,y], [u,v], [[0,1], [1,-1]], [y=3]);

bproc1 := proc(d) 
	if d = 1 then
		return 3*u+ v;
	else 
		return 0;
	end if; 
end proc;  

pso := PowerSeries(bproc1, variables={u,v});
puso9 := PuiseuxSeries(pso, [x,y], [u,v], [[0,1], [1,-1]], [y=3]);

bproc2 := proc(d) 
	if d = 1 then
		return 3*u+ v;
	elif d = 2 then
		return u*v + u^2;
	elif d mod 2 = 0 then
		return 0;
	else  
		return u^(d-1)*v;
	end if; 
end proc;  

pso := PowerSeries(bproc2, variables={u,v});
puso10 := PuiseuxSeries(pso, [x,y], [u,v], [[0,1], [1,-1]], [x=-1,y=3]);

# Test converto.
Try[testnoerror]("test conv 1", puso1:-ConvertToPSO(puso1), 'assign'='pso1');
Try[testerror]("test conv 2", puso2:-ConvertToPSO(puso2));
Try[testnoerror]("test conv 3", puso3:-ConvertToPSO(puso3), 'assign'='pso3');
Try[testnoerror]("test conv 4", puso4:-ConvertToPSO(puso4), 'assign'='pso4');
Try[testnoerror]("test conv 5", puso5:-ConvertToPSO(puso5), 'assign'='pso5');
Try[testnoerror]("test conv 6", puso6:-ConvertToPSO(puso6), 'assign'='pso6');
Try[testerror]("test conv 7", puso7:-ConvertToPSO(puso7), 'assign'='pso7');
Try[testnoerror]("test conv 8", puso8:-ConvertToPSO(puso8), 'assign'='pso8');
Try[testnoerror]("test conv 9", puso9:-ConvertToPSO(puso9), 'assign'='pso9');
Try[testnoerror]("test conv 10", puso10:-ConvertToPSO(puso10), 'assign'='pso10');

# Test Truncate.
poly := Truncate(puso1, 10, mode=absolute);
Try[verify,normal]("test trunc 1", pso1:-Truncate(pso1, 10), poly);

poly := Truncate(puso3, 10, mode=absolute);
Try[verify,normal]("test trunc 3", pso3:-Truncate(pso3, 10), poly);

poly := Truncate(puso4, 10, mode=absolute);
Try[verify,normal]("test trunc 4", pso4:-Truncate(pso4, 10), poly);

poly := Truncate(puso5, 10, mode=absolute);
Try[verify,normal]("test trunc 5", pso5:-Truncate(pso5, 10), poly);

poly := Truncate(puso6, 10, mode=absolute);
Try[verify,normal]("test trunc 6", pso6:-Truncate(pso6, 10), poly);

poly := Truncate(puso8, 10, mode=absolute);
Try[verify,normal]("test trunc 7", pso8:-Truncate(pso8, 10), poly);

poly := Truncate(puso9, 10, mode=absolute);
Try[verify,normal]("test trunc 8", pso9:-Truncate(pso9, 10), poly);

Try[testerror]("test trunc ", pso10:-Truncate(pso10, 10));

#end test