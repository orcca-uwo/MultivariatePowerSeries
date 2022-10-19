#test
with(TestTools):
Try[testnoerror]("test 0", with(MultivariatePowerSeries));

random_seed := randomize();
randomize(random_seed);
PrintOnFail("seed", random_seed);

# We generate a random zero ray of dimension n with  
# entries from equal to p to q. Note, the ray could not be grevlex positive.
RandomZeroRay := proc(n :: nonnegint, 
                                   p :: negint, 
                                   q :: nonnegint,$)
    local r := Array(1 .. n);
    local s := 0;
    local roll := rand(p..q);

    for local i  from n to 2 by -1 do 
        r[i] := roll();
        s :=  s + r[i];
    end do;

    r[1] := -s;

    return convert(r, ':-list');
end proc;

# We generate a random ray of dimension n and wight of w>0.  
RandomRay := proc(n :: nonnegint, 
                                    w :: posint, $)
    local r := Array(1 .. n);
    local s := 0;
    local roll := rand(-w..w);

    for local i  from n to 2 by -1 do 
        r[i] := roll();
        s :=  s + r[i];
    end do;

    r[1] := -s+w;

    return convert(r, ':-list');
end proc;

kernelopts(opaquemodules=false):
PuSO := MultivariatePowerSeries:-PuiseuxSeriesObject:
MakeRaysCompatible := MultivariatePowerSeries:-PuiseuxSeriesObject:-MakeRaysCompatible:
kernelopts(opaquemodules=true):

ord := [x1, x2, x3, x4, x5];
pso := PowerSeries(1/(1+u1*u2*u3*u4));
e := [x1=-5, x2=2];
rays := [[-2,1,2,-1,0], [3,2,1,2,3], [2,-2,-4,5,-1], [-3,4,0,-1,0]];
ordCV := [u1, u2, u3, u4];

Try[testnoerror]("test 1", PuSO(pso, ord, ordCV, rays, e), 'assign'='puso4');
newOrd := puso4:-GetPuiseuxSeriesOrder(puso4);

#### zero-ray test.
for j from 1 to 100 do
	R := [];
	for i from 1 to 8 do
		r := RandomZeroRay(5, -4, 10);
		if PuSO:-Positive(PuSO, r) then
			R := [op(R), r];
		end if;
	end do;

	if R <> [] then
		Try[testnoerror](cat("test ", j, " 2"), MakeRaysCompatible(puso4, R, newOrd), 'assign'='newRays');

		M := Matrix(newRays);
		raysAsVectors := map(convert, R, ':-Vector');
		result := [seq(LinearAlgebra:-LinearSolve(M^+, r), r in raysAsVectors)];
		Try(cat("test ", j, " 3"), andseq(PuSO:-Positive(puso4, r), r in result), true);

		Try(cat("test ", j, " 4"), andseq(type(convert(r,list),list(nonnegint)), r in result), true);
	end if;
end do;

### positive weight rays.
roll := rand(1..10);
for j from 1 to 100 do
	R := [];
	for i from 1 to 8 do
		r := RandomRay(5, roll());
		R := [op(R), r];
	end do;

	if R <> [] then
		Try[testnoerror](cat("test ", j, " 5"), MakeRaysCompatible(puso4, R, newOrd), 'assign'='newRays');

		M := Matrix(newRays);
		raysAsVectors := map(convert, R, ':-Vector');
		result := [seq(LinearAlgebra:-LinearSolve(M^+, r), r in raysAsVectors)];
		Try(cat("test ", j, " 6"), andseq(PuSO:-Positive(puso4, r), r in result), true);

		Try(cat("test ", j, " 7"), andseq(type(convert(r,list),list(nonnegint)), r in result), true);
	end if;
end do;
#end test