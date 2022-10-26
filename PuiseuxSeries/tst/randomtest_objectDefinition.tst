#test 200
with(TestTools):
Try[testnoerror]("test 0", with(MultivariatePowerSeries));

random_seed := randomize();
randomize(random_seed);
PrintOnFail("seed", random_seed);

kernelopts(opaquemodules=false):
PuSO := MultivariatePowerSeries:-PuiseuxSeriesObject:
MakeRaysCompatible := MultivariatePowerSeries:-PuiseuxSeriesObject:-MakeRaysCompatible:
kernelopts(opaquemodules=true):

###########################
## CreatePuSOFromRatFunc ##
###########################
mp := [u=x^(1/2), v=x^(1/3)*y^(-1/6), w=y^(1/2)*z^(-1/2)]; 
roll := rand(-10..10);
for i from 1 to 10 do
	num := randpoly([u,v,w],degree=4);
	den := randpoly([u,v,w],dense,degree=4);
	e := [x=roll(), y=roll()];

	Try[testnoerror]("test ratfun 1", PuSO(num/den, mp, e), 'assign'='pus_rat');

	Try[testnoerror]("test ratfun 2", PuSO:-Inverse(pus_rat), 'assign'='inv');
	Try[testnoerror]("test ratfun 3", PuSO:-BinaryMultiply(pus_rat, inv), 'assign'='pusMult');
	Try("test ratfun 4", PuSO:-Truncate(pusMult, 40), 1);
end do;

#end test