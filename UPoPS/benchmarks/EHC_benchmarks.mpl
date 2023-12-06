# BenchMarks EHC.
with(MultivariatePowerSeries);


my_test := proc(p, v, vars, D, print_output, multivariate:=false,
				returnleadingcoefficient:= false)
	local u := UnivariatePolynomialOverPowerSeries(p, v);

	## New EHC.
	local t1 := 0;
	local EHC1;
	t1, EHC1 := CodeTools:-Usage(u:-ExtendedHenselConstruction(u, a)
	                              , output = ['cputime', 'output']);

	local exp, pol;
	local t := 0;
	t, exp := CodeTools:-Usage([seq(Truncate(pol, D), pol in EHC1)]
	                              , output = ['cputime', 'output']);
	t1	 += t;
	if print_output then
		print("New EHC ", exp);
	end if;

	## Checking the iteration reached in EHC.
	local c;
    if returnleadingcoefficient then
        c := GetCoefficient(EHC1[2],Degree(EHC1[2]));
    else
        c := GetCoefficient(EHC1[1],Degree(EHC1[1]));
    end if;
    local c_pso := c:-GetPowerSeries(c);
    local state := c_pso:-GetAncestor(c_pso, ':-state');
    local k := state:-currentJ;

	## Puiseux factorization.
	local t2 := 0;
	local PF;
	t := 0;
	if multivariate=false then
		t2, PF := CodeTools:-Usage(PuiseuxFactorize(u)
		                              , output = ['cputime', 'output']);

	    t, exp := CodeTools:-Usage([seq(Truncate(pol, D), pol in PF)]
	                       	       , output = ['cputime', 'output']);
	end if;
	t2 += t;
	if print_output then
		print("Puiseux theorem ", exp);
	end if;

	##
	# Calling the RegularChains code.
	local U := RegularChains:-UnivariatePolynomialOverPowerSeries(vars, v):;
	local pnt := 0;
	if multivariate=true then
		pnt := [0,0];
	end if;
	local iter := k;
	local OutputFlag::name:='parametric';
	local parametricVar::name:= T;
	local verificationFlag::boolean:=false;

	local t3 := 0;
	local EHC2;
	t3, EHC2 :=  CodeTools:-Usage(U:-ExtendedHenselConstruction(p, pnt, iter, OutputFlag, parametricVar, verificationFlag)
	                              , output = ['cputime', 'output']);

	if print_output then
		print("Old EHC ", EHC2);
	end if; 

	return [t1, t2, t3];
end proc:

# Take into account that We create the new EHC and Puiseux factorization 
# each time that we call my_test. So, in reality, we are even faster.

print_output := false;

### 1 polynomial.
p := x^5 + x^4*y - 2*x^3*y - 2*x^2*y^2 + x*(y^2 - y^3) + y^3;
ts := my_test(p, x, [y], 5, print_output); 

print("k", 5, "t1: ", ts[1], "t2: ", ts[2], "t3: ", ts[3]);

ts := my_test(p, x, [y], 15, print_output); 

print("k", 15, "t1: ", ts[1], "t2: ", ts[2], "t3: ", ts[3]);

ts := my_test(p, x, [y], 30, print_output); 

print("k", 30, "t1: ", ts[1], "t2: ", ts[2], "t3: ", ts[3]);

### 2 polynomial.
p := x^3 -x^2*y^2 -x*y^3 + y^4;
ts := my_test(p, x, [y], 5, print_output); 

print("k", 5, "t1: ", ts[1], "t2: ", ts[2], "t3: ", ts[3]);

ts := my_test(p, x, [y], 15, print_output); 

print("k", 15, "t1: ", ts[1], "t2: ", ts[2], "t3: ", ts[3]);

ts := my_test(p, x, [y], 30, print_output); 

print("k", 30, "t1: ", ts[1], "t2: ", ts[2], "t3: ", ts[3]);

### 3 polynomial.
p := y^3*(y+1)*x^2 + RootOf(x^2+1)*x*y + y^2;
ts := my_test(p, x, [y], 5, print_output, false, true); 

print("k", 5, "t1: ", ts[1], "t2: ", ts[2], "t3: ", ts[3]);

ts := my_test(p, x, [y], 15, print_output, false, true); 

print("k", 15, "t1: ", ts[1], "t2: ", ts[2], "t3: ", ts[3]);

ts := my_test(p, x, [y], 30, print_output, false, true); 

print("k", 30, "t1: ", ts[1], "t2: ", ts[2], "t3: ", ts[3]);

### 4 polynomial.
p := x^4 + 3*(y+y^2)*x^3 + 3*(y^2+y^3)*x^2+y^3*x;
ts := my_test(p, x, [y], 5, print_output); 

print("k", 5, "t1: ", ts[1], "t2: ", ts[2], "t3: ", ts[3]);

ts := my_test(p, x, [y], 15, print_output); 

print("k", 15, "t1: ", ts[1], "t2: ", ts[2], "t3: ", ts[3]);

ts := my_test(p, x, [y], 30, print_output); 

print("k", 30, "t1: ", ts[1], "t2: ", ts[2], "t3: ", ts[3]);


### 5 polynomial.
p := (x-1)*(x-2)*(x-3)+ y*(x^2+x);
ts := my_test(p, x, [y], 5, print_output); 

print("k", 5, "t1: ", ts[1], "t2: ", ts[2], "t3: ", ts[3]);

ts := my_test(p, x, [y], 15, print_output); 

print("k", 15, "t1: ", ts[1], "t2: ", ts[2], "t3: ", ts[3]);

ts := my_test(p, x, [y], 30, print_output); 

print("k", 30, "t1: ", ts[1], "t2: ", ts[2], "t3: ", ts[3]);

### 6 polynomial.
p := (x-1)*(x-2)*(x-3)*(x-4)+ y*(x^3+x);
ts := my_test(p, x, [y], 5, print_output); 

print("k", 5, "t1: ", ts[1], "t2: ", ts[2], "t3: ", ts[3]);

ts := my_test(p, x, [y], 15, print_output); 

print("k", 15, "t1: ", ts[1], "t2: ", ts[2], "t3: ", ts[3]);

ts := my_test(p, x, [y], 30, print_output); 

print("k", 30, "t1: ", ts[1], "t2: ", ts[2], "t3: ", ts[3]);

### 7 polynomial.
p := x^3 + (y - z + z^2)*x^2 - (y + z + y^2 - z^2)*x + (y^2 - z^3);
ts := my_test(p, x, [y, z], 5, print_output, true); 

print("k", 5, "t1: ", ts[1], "t2: ", ts[2], "t3: ", ts[3]);

ts := my_test(p, x, [y, z], 15, print_output, true); 

print("k", 15, "t1: ", ts[1], "t2: ", ts[2], "t3: ", ts[3]);
