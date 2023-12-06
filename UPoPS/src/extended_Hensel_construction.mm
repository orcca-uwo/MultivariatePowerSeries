# To compute the first non-zero term of a PSO up to 
# bound.
local firstNonzeroTerm::static := proc(s::PowerSeriesObject, 
										bound::integer, $)
	local ana := s:-GetAnalyticExpression(s);

    # We check if we have a polynom or not.
	if type(ana, ':-polynom') and ana <> undefined then
		local as_sum := convert(ana, ':-list', :-`+`);
        # We compute the degree of all the monomials.
		local degs := map(degree, as_sum);

		return min(degs);
	else 
        # We try to find a non-zero homogeneous part
        # up to degree bound.
        for local i from 0 to bound do
          if s:-HomogeneousPart(s, i) <> 0 then
            return i;
          end if;
        end do;

        # Did not find a nonzero value before the bound.
        return FAIL;
	end if;

end proc;

# Given a monomial, a list of 
# variables, and a main variable, returns a list with 
# two entries: the degree in the main variable and the 
# (total) degree in the other variables.
local exponents::static := proc(expr::polynom, vars::list(name), $)
    local result := map2(degree, expr, vars);

    return [result[1], add(result[2..])];
end proc;

# For a given polynomial we compute a list of its monomial called as_sum.
# Then, for each monomial in as_sum, we compute 
# a list of pairs [exponent, sum_exponent], where 'exponent' is the
# main degree and sum_exponent is the total degree - 'exponent'. 
local exponentList::static := proc(self::UnivariatePolynomialOverPowerSeriesObject,
									poly::polynom, $) 
	local poly_normal := normal(poly);
	local as_sum := convert(poly_normal, ':-list', :-`+`);
	local mvar := self:-MainVariable(self);
	local vars := ListTools:-MakeUnique([mvar, op(self:-Variables(self))]);
	
	return map(self:-exponents, as_sum, vars);
end proc;

# Returns the slope between two points point1=(x_1,y_1 )
# and point2= (x_2,y_2 ) in the plane.
local slope::static := proc(point1::[integer, integer],
							point2::[integer, integer], $)
    return((point1[2] -  point2[2] )/(point1[1] -  point2[1]));
end proc:

# Computes the newton line of a monic square-free polynomial
# self.
local NewtonLine::static := proc(self :: UnivariatePolynomialOverPowerSeriesObject,
                               Hensel_bound::integer, $)
    local d := DEGREE(self);
    local v := self:-MainVariable(self);
    local poly := self:-GetAnalyticExpression(self);
    local term, slopes, min_slope, delta, e;
	local d_as_point := [d, 0];

    if type(poly, ':-polynom') and poly <> undefined then
        # If self has a polynomial analytic expression,
        # then we compute a list of all the exponents of
        # the monomials of self.
    	local exp_ls := self:-exponentList(self, poly);

        # Then we compute all the possible slops of the
        # lines from the points in exp_ls to d_as_point.
    	slopes:= [seq(ifelse(e<>d_as_point,
    					self:-slope(d_as_point, e), NULL), e in exp_ls)];

        # We compute the newton slope.
        min_slope:= max(slopes);

    else 
        # We compute the first non-zero term in upoly[0] 
        # up to bound and compute its degree. 
    	local upoly := self:-upoly;
    	local sm_deg := self:-firstNonzeroTerm(upoly[0], Hensel_bound);

        # If firstNonzeroTerm did not fail, then
        # min_slope is computed.
    	if sm_deg <> FAIL then 
    		min_slope := self:-slope(d_as_point, [0, sm_deg]);
    	else 
    		min_slope := -infinity;
    	end if;

    	local new_slope := -infinity;
    	for local i from 1 to d-1 do 
    		local my_bound := min(Hensel_bound, floor((i - d) * min_slope));
    		
            # We compute the first non-zero term in upoly[i] 
            # up to bound and compute its degree. 
            sm_deg := self:-firstNonzeroTerm(upoly[i], my_bound);
    		
            # If firstNonzeroTerm did not fail, then
            # this means we have a point below the candidate
            # to newton line. Hence, we update new_slope.
            if sm_deg <> FAIL then 
	    		new_slope := self:-slope(d_as_point, [i, sm_deg]);
	    	else 
	    		new_slope := -infinity;
	    	end if;

	    	min_slope:= max([min_slope, new_slope]);

    		# In case we reach the minimum possible value.
            if min_slope = 0 then
    			break;
    		end if;
    	end do;

    end if;
    delta:= min_slope*(-d) +0;

    return d, delta;
end proc;


# Returns the newton polynomial of self together 
# with d and delta.
local NewtonPolynomial::static:= proc(self::UnivariatePolynomialOverPowerSeriesObject, 
										Hensel_bound::integer, 
                                        v::name, u::name, {multivariate :: truefalse := false}, $)
    local d, delta;
	# First, we compute the newton line of self.
	# Note d is integer and delta rational.
	d, delta := self:-NewtonLine(self, Hensel_bound);

    if delta=infinity then 
        error "it is not possible to compute the Newton line with the current bound %1, try to increase it", Hensel_bound;
    end if;

    local slope := normal(delta/d);
    
    # delta*(degree(MainVarList[i],v) )+d*( degree(coeffList[i][j],u) ) = delta*d
	local terms := 0;
    local mono, vars_as_conts;
    if multivariate then 
        local var;
        local vars := self:-Variables(self) minus {v};
        vars_as_conts := [seq(var=var(0), var in vars)];
    end if; 

    for local i from 0 to d do 
        local deg := (d - i) \*slope;
        if type(deg, ':-integer') then # remove this condition once we deal with Puiseux series!
            mono := PowerSeriesObject:-Truncate(self:-upoly[i], deg);
            if multivariate then 
                mono := subs(vars_as_conts, mono)*u^deg;
            end if;
            terms += mono * v^i;
        end if;
    end do;

    local newton_pol := add(terms);

	return d, delta, newton_pol;
end proc;

local my_factor::static := proc(p :: depends(polynom(anything, x)),
                                 x :: name, $)
  if hastype(p, ':-And'(':-constant', ':-Not'(':-radalgnum'))) then
    local sols := solve(p, [x]);
    ASSERT(type(sols, ':-list'([':-identical'(x) = ':-anything'])));
    local facs := map(lhs - rhs, map(op, sols));
    facs := map(simplify, facs);
    local f, z;
    return [seq([f, add(ifelse(z = f, 1, 0), z in facs)], f in {op(facs)})];
  else
    return evala(PolynomialTools:-Splits(p, x)[2]);
  end if;
end proc;

# Create the list of initial factors for extended Hensel construction.
local initialFactors::static := proc(self::UnivariatePolynomialOverPowerSeriesObject,
										Hensel_bound::integer, 
                                        v::name, u::name, {multivariate :: truefalse := false}, $)
	local d, delta, newton_pol;
	# We compute the newton polynomial.
    d, delta, newton_pol := self:-NewtonPolynomial(self, Hensel_bound, v, u, _options['multivariate']);
	#local v := self:-MainVariable(self);
	#local u := mul((self:-Variables(self) minus {v}));

	# Computing d_hat and delta_hat.
	local fract := delta/d;
    local d_hat:= denom(fract);
    local delta_hat:= numer(fract);

	# We evaluate the variable u=1. Now poly is univariate.
    local poly:= eval(newton_pol, u=1);

    # Univariate complete factorization of poly.
    # NOTE: We could consider just calling 
    # PolynomialTools:-Split and then pairwise 'is' on 
    # any factors that have non-algebraic constants in them.
    local Facs:= my_factor(poly ,v);

    return d_hat, delta_hat, Facs;
end proc;

# Computes the extended gcd among the polynomials in polys
# and returns the my_gcd and 'Bezou_coeff', where  my_gcd is the gcd(polys)
# and 'Bezou_coeff' are the bezout coefficient w.r.t polynomials in polys.
local extendedGcd::static := proc(polys::list(polynom), v::name, $)
	local size := numelems(polys);
	ASSERT(size>=2, "first parameter must has at least two polynomials");
	
	local Bezou_coeff := Array(1..size, fill=1);

	local my_gcd := polys[1];
	local quo1, quo2, dummy, p1, p2;

    for local i from 2 to size do
        p1 := subsindets(my_gcd, 'anything'^'fraction', convert, RootOf);
        p2 := subsindets(polys[i], 'anything'^'fraction', convert, RootOf);
        my_gcd := evala(Gcdex(p1, p2, v, 'quo1', 'quo2'));

        for local j from 1 to i-1 do
            Bezou_coeff[j]:=  AUTO_EXPAND(dummy, quo1*Bezou_coeff[j]);
        end do;
        Bezou_coeff[i]:= AUTO_EXPAND(dummy, quo2*Bezou_coeff[i]);
    end do;

    return my_gcd, Bezou_coeff;
end proc:

# Returns a array containing Lagrange Interpolation
# Polynomials w.r.t 'polys_in'.
# Note that we use the univariate polynomials polys_in. In
# the original algorithm the newton polynomial is lifted,
# then the second variable is evaluated in 1 again. In this 
# version, we do not lift the newton polynomial.
local LagrangeInterpolationPolynomials::static := proc(polys_in::list([algebraic, nonnegint]), 
												v::name, u::name, delta_hat::nonnegint, 
                                                multivariate::truefalse := false, $)
	# polys_in is a list of univariate polynomials in v.
	local size := numelems(polys_in);
	
	# We compute the mul(polys_in)/polys_in[i] for i=1..size.
	# Note that first we expand polys_in[i][1]^polys_in[i][2].
	local i, j, dummy;
	local polys := [seq(polys_in[i][1]^polys_in[i][2], i=1..size)];
    local poly_list:= [seq(AUTO_EXPAND(dummy, mul(polys[j], j=1..i-1 )* mul(polys[j], j=i+1..size )),
                                         i=1..size  )];

    # Computing the extended gcd of polynomials in poly_list.
    # We know for sure the gcd must be 1.
    local my_gcd, Bezou_coeff;
    
    my_gcd, Bezou_coeff:= extendedGcd(poly_list, v);

	ASSERT(my_gcd=1, "gcd must be 1");

    # d = deg_x(mul(polys_in)).
    # Note that we use the univariate polynomials polys_in.
    # There is no need to lift these polynomials yet.
	local f;
    local d:= add(degree(f[1], v)*f[2], f in polys_in);

	# The output.
	local LagIntPoly:= Array(0..(d-1));

    # Here we aim at computing the Lagrange Interpolation Polynomials.
    for local l from 0 to d-1 do
        local QuoTable:= Array(1..size);
        local RemTable:= Array(1..size);
        
        # We compute Qi and Ri, such that
        # Bezou_coeff[i]*v^l = Qi*polys[i] + Ri 
        # for l=0..d-1. 
        for i from 1 to size - 1 do
            # Here we check whether
            # degree(v^l*BezouTable[i])>degree(polys[i]) or not.
            if degree(Bezou_coeff[i],v) + l >= degree(polys[i],v) then
                local q;

                RemTable[i]:= evala(Rem(Bezou_coeff[i]*v^l, 
                                                    polys[i],v,'q'));
                QuoTable[i]:= q;
            else
                RemTable[i]:= Bezou_coeff[i]*v^l;
                QuoTable[i]:= 0;
            end if;
        end do;

        # We lift our polynomials back.
        local last_term := Groebner:-Homogenize(
                        AUTO_EXPAND(dummy, Bezou_coeff[size]*v^l + add(QuoTable[i], i=1..size-1)*polys[size]),u);
        local temp:= [seq(Groebner:-Homogenize(RemTable[i],u), i=1..size-1), 
        			  last_term];

        temp:= [seq(AUTO_EXPAND(dummy, temp[i]*u^(degree(polys[i],v)-degree(temp[i],v) )),
            			i=1..size)];

        LagIntPoly[l] := map(eval, temp, u= u^(delta_hat));
    end do;

    ASSERT(andseq(andseq(degree(LagIntPoly[j][i], v) < degree(polys[i], v), i=1..size), 
    				j=0..d-1),
    		"degree in main variable of the i-Lagrange polynomial must be smaller than"
    		" the degree in the main variable of G_i");
    
    ASSERT(andseq(andseq(PolynomialTools:-IsHomogeneous(LagIntPoly[j][i]), i=1..size), 
    				j=0..d-1),
    				"Lagrange polynomials must be homogeneous");

    return LagIntPoly;
end proc:

# Returns a array containing Lagrange Interpolation
# Polynomials w.r.t 'polys_in'.
# Note that we use the univariate polynomials polys_in. In
# the original algorithm the newton polynomial is lifted,
# then the second variable is evaluated in 1 again. In this 
# version, we do not lift the newton polynomial.
local FastLagrangeInterpolationPolynomials::static := proc(polys_in::list([algebraic, nonnegint]), 
                                                v::name, u::name, delta_hat::nonnegint, $)
    # polys_in is a list of univariate polynomials.
    local size := numelems(polys_in);

    # We compute the mul(polys_in)/polys_in[i] for i=1..size.
    # Note that first we expand polys_in[i][1]^polys_in[i][2].
    local i, j, l, k;
    local polys := [seq(polys_in[i][1]^polys_in[i][2], i=1..size)];
    local dummy;
    local poly_list:= [seq(AUTO_EXPAND(dummy, mul(polys[j], j=1..i-1)* mul(polys[j], j=i+1..size )), 
                                        i=1..size  )];

    # Note that we use the univariate polynomials polys_in.
    # There is no need to lift these polynomials yet.
    local f;
    local deg:= add(degree(f[1], v)*f[2], f in polys_in);

    # The output.
    local LagIntPoly :=  Array(0..deg-1, [seq(Array(1..size), i=0..deg-1)]);

    for i from 1 to size do
        # The idea is to solve a system of matrix 
        # equations. We want to compute X such that
        # W*X=B, where W is the Wronskian Matrix and B
        # the vector of derivatives from the theorem.
        # We also know that the inverse of W is the
        # product of the matrices M_1 and M_2 computed
        # bellow. 
        local mult:= polys_in[i][2];
        local vv:= -(polys_in[i][1] - v);

        # We compute the M_1_j's matrices.
        local M_1_j:= [seq(Matrix(
                        [seq([seq(0, k=0..j-1), 1, 
                            seq(0, k=j+1..mult-1)], j=0..l-1), 
                        [seq(0, k=0..l-1), 1/(factorial(l)*poly_list[i]), 
                         seq(0, k=l+1..mult-1)], 
                         seq([seq(0, k=0..l-1), (-1)*binomial(j,l)*(diff(poly_list[i],v$(j-l)))/(poly_list[i]),  
                         seq(0, k=l+1..j-1), 1, seq(0, k=j+1..mult-1)], j=l+1..mult-1)]), l=0..mult-1)];

        # Since poly_list are univariate, we just need to evaluate in v=vv.
        local M_1_eval:= eval(M_1_j, v=vv);

        # Compute the multiplication of all matrices in the list M_1_eval.
        local M_1:= M_1_eval[-1];
        for j from mult-1 to 1 by -1 do
            LinearAlgebra:-Multiply(M_1, M_1_eval[j], inplace);
        end do;

        # Compute the M_2 matrix.
        # Note: In the old version, they were evaluating M_2 in u=1.
        # This is clearly not necessary. 
        local M_2:= Matrix(eval([seq([seq(0, k=0..j-1),
                                seq((-1)^(j+k)*binomial(k,k-j)*v^(k-j), k=j..mult-1)], j=0..mult-1)], v= vv));

        # Compute the inverse of Wronskian Matrix. 
        # The theory says the inverse of Wronskian Matrix is
        # the product of M_1 and M_2.
        local inv_wrons:= LinearAlgebra:-Multiply(M_2, M_1);

        # For each (l) compute the Lagrange polynomial W_i^(l).
        for  l from 0 to deg-1 do
            local B_matrix:= Matrix(eval(eval([[v^(l)*u^(deg-l)],
                                seq([diff(v^(l)*u^(deg-l), v$j)], j=1..mult-1)], v=vv), u=1));
            local answer:= LinearAlgebra:-Multiply(inv_wrons, B_matrix);

            local my_expd:= AUTO_EXPAND(dummy, add(answer[j][1]*v^(j-1)*u^(mult-j+1), j=1..mult));
            LagIntPoly[l][i] :=  map(eval, my_expd, u= u^(delta_hat));
        end do;
    end do;

    return LagIntPoly;
end proc;

# Computes the required coefficients for doing one lifting
# in Hensel construction method.
# Remark:  All the terms
# of 'poly' are of the form 'v^(d-i) * u^((k+i*delta_hat)/d_hat)'
# where i=0..d and 'k' is fixed(it is passed as an argument to
# function ComputeCoeffs). Then for each 'i' it returns the c_i^(k)*u^k
# where  c_i^(k) is the coefficient of the term 'v^(d-i)*u^((k+i*delta_hat)/d_hat)'.
local computeCoeffs::static := proc(self::UnivariatePolynomialOverPowerSeriesObject,
                                    poly::polynom, 
                                    v::name,
                                    u::name,
                                    d::integer, 
                                    d_hat::integer, 
                					delta_hat::integer, 
                                    k :: integer, $)
    local k2:= k+1;

    # Initializing the coefficient array.
    local CoeffTable:= Array(0..d-1);

    for local l from 0 to d-1 do
        local temp:= Groebner:-TrailingTerm(coeff(poly, v, l), plex(v, u));

        # Anything with degree greater or equal than
        # k2 +(d-l)*delta_hat must be zero in the current
        # iteration, currentJ, of the Hensel algorithm.
        if degree(temp[2], u) >= k2 +(d-l)*delta_hat then
            CoeffTable[l]:= 0;
        else
            CoeffTable[l]:= temp[1];
        end if;
    end do;

    return CoeffTable;
end proc:

# To be computed: The product of the iteration
# k+1=currentJ.
# Paper formulas:
# Assume Y=u_hat=u^(1/d_hat).
### Equality (1) ###
# DELTA_j^k*Y^k = DELTA(G)_j^k.
# DELTA_j^0 = G_j^0.
### Equality (2) ###
# P_{j-1}^k = p_{j-1}^{k,0} + p_{j-1}^{k,1}*Y + ... + p_{j-1}^{k,k}*Y^k .
### k+1 = 1 ###
# P_j^1 := P_{j-1}^1*G_j^0, P_1^1 := G_1^0, j>0.
### k+1 > 1 ###
## P_2 formulas ##
# q_2^{k+1}*Y^k := (DELTA_1^0*DELTA_2^k + DELTA_1^0*DELTA_2^k)*Y^k.
# P_2^{k+1} := P_2^k + q_2^{k+1}*Y^k + sum(DELTA_1^i*DELTA_2^{k+1-i}*Y^{k+1}, i=1..k).
# Thanks to equality (1)
# q_2^{k+1}*Y^k := DELTA(G)_1^0*DELTA(G)_2^k + DELTA(G)_1^0*DELTA(G)_2^k.
# P_2^{k+1} := P_2^k + q_2^{k+1}*Y^k + sum(DELTA(G)_1^i*DELTA(G)_2^{k+1-i}, i=1..k).
## P_j formulas, j>2 ##
# q_j^{k+1}*Y^k := (p_{j-1}^{k+1,0}*DELTA_j^k + q_{j-1}^{k+1}*DELTA_j^0)*Y^k.
# P_j^{k+1} := P_j^k + q_j^{k+1}*Y^k + sum(p_{j-1}^{k+1,i}*DELTA_j^{k+1-i}*Y^{k+1}, i=1..k+1).
# Thanks to equality (1)
# q_j^{k+1}*Y^k := p_{j-1}^{k+1,0}*DELTA(G)_j^k + (q_{j-1}*Y^k)^{k+1}*DELTA_j^0.
# P_j^{k+1} := P_j^k + q_j^{k+1}*Y^k + sum(p_{j-1}^{k+1,i}*Y^i*DELTA(G)_j^{k+1-i}, i=1..k+1).
# Note: fast_product[j] stores the terms
# of P[j] at the k iteration, ie P[j][i]=p_j^{k,i}*Y^i, i=0,...,k.
# k starts at 0. 
# At each iteration, we modify P[j][k] and we also add the
# next component P[j][k+1]. 
# fast_product[size+1] stores P_size^{k+1}.
local FastProduct::static := proc(state::record, $)
    local k := state:-currentJ-1;
    local fast_product := state:-fast_product;
    local upolys := state:-upolys;
    local multivariate := state:-multivariate;

    local init_upop := state:-init_upop;
    local size := numelems(upolys);
    if state:-lc_included then 
        size -=1;
    end;

    # This is equal to to G_i^(0),
    # with i=1..size.
    local LiftedFactors := state:-LiftedFactors;
    local d_hat := state:-d_hat;
   
    # Base case.
    local i, dummy;
    if k+1=1 then 
        fast_product[1] ,= LiftedFactors[1][0], 0;
        for local j from 2 to size do
            fast_product[j] ,= fast_product[j-1][0]*LiftedFactors[j][0], 0;
        end do;

        fast_product[size+1] := AUTO_EXPAND(dummy, fast_product[size][0]);

        return fast_product[size+1];
    end if;

    # First factor.
    # G_1^k.
    fast_product[1] ,= LiftedFactors[1][k];

    # Second factor.
    local my_sum := add(LiftedFactors[1][i]*LiftedFactors[2][k+1-i], i=1..k);
    local recursive_q := LiftedFactors[1][0]*LiftedFactors[2][k]
                                + LiftedFactors[1][k]*LiftedFactors[2][0];

    fast_product[2][k] := fast_product[2][k] + recursive_q;
    fast_product[2] ,= my_sum;

    for local j from 3 to size do 
        recursive_q:= fast_product[j-1][0]*LiftedFactors[j][k] 
                        + recursive_q*LiftedFactors[j][0];

        my_sum := add(fast_product[j-1][i]*LiftedFactors[j][k+1-i], i=1..k+1);
        
        fast_product[j][k] := fast_product[j][k] + recursive_q;
        fast_product[j] ,= my_sum;
    end do;
    fast_product[size+1]:= AUTO_EXPAND(dummy, fast_product[size+1]+recursive_q+fast_product[size][k+1]);


    return fast_product[size+1];
end proc:


# It runs 1 iteration of the Hensel algorithm.
# It takes care of increasing the approximation of the
# Hensel factors (upolys).
local MainIteration ::static := proc(state::record, $)
    local init_upop := state:-init_upop;
    local d := DEGREE(init_upop);
    local upolys := state:-upolys;
    # Iteration to be run.
    local currentJ := state:-currentJ;

    local LiftedFactors := state:-LiftedFactors;

    local multivariate := state:-multivariate;
    local v := state:-v;
    local u := state:-u;

    # Product of the current factorization.
    local fac_prod := 1;
    # Note: We do not want to include the upolys[0] (the lc
    # of the input polynomial).
    local upolys_size := numelems(upolys);
    if state:-lc_included then 
        upolys_size -=1;
    end;

    # Fast product of the factors up to the iteration CurrentJ-1.
    if upolys_size>1 then
        fac_prod := init_upop:-FastProduct(state);
    else 
        fac_prod := upolys[1]:-Truncate(upolys[1]);
    end if; 

    # Initial polynomial.
    local d_hat := state:-d_hat;
    local delta_hat := state:-delta_hat;

    # We get all the non-zero terms of init_upop mod S_{k+1}.
    local id;
    local init_upop_poly := add(function:-Truncate(init_upop:-upoly[id],
                                                     currentJ+ (d-id)*delta_hat)*x^id, id=0..d);

    if multivariate then 
        local var;
        local vars := init_upop:-Variables(init_upop) minus {v};
        local vars_as_conts := [seq(var=var(0)*u, var in vars)];
        init_upop_poly := subs(vars_as_conts, init_upop_poly);
    end if;
    init_upop_poly := eval(init_upop_poly, u=u^d_hat);

    # The error.
    local dummy;
    # The algorithm does not work without this simplification (Algebraic:-Normal).
    local delta_F := Algebraic:-Normal(init_upop_poly - fac_prod);

    # Compute the li's.
    # All the output must be multiplied by u^(currentJ+1).
    local l_i_s := init_upop:-computeCoeffs(init_upop, delta_F, v, u, d, d_hat, 
                                            delta_hat, currentJ);

    l_i_s := l_i_s*u^(currentJ);
    # Multiply the Lagrange polys by the li's.
    local lag_polys := state:-lag_polys;

    # The delta_G_i's of the theorem.
    local poly_increase := Array(1..upolys_size, 0);
    for local i from 0 to d-1 do 
        for local j from 1 to upolys_size do 
            poly_increase[j] := poly_increase[j] + 
                                lag_polys[i][j]*l_i_s[i]; 
        end do;
    end do;


    for local j from 1 to upolys_size do 
        # it does not work in the multivariate case
        poly_increase[j] := AUTO_EXPAND(dummy, poly_increase[j]); 
        ASSERT(degree(poly_increase[j], u)<=d*delta_hat+currentJ);

        LiftedFactors[j] ,= poly_increase[j];
    end do;
   
    # Update all the factors.
    # Note that we also update the coefficients that
    # must be zero.
    for local i from 1 to upolys_size do 
        local m_i := DEGREE(upolys[i]);
        local coefficients := PolynomialTools:-CoefficientList(poly_increase[i], v);
        local size := numelems(coefficients);
        
        for local j from 0 to m_i-1 do
            
            local coefficient := 0;
            if numelems(coefficients)>=j+1 then
                coefficient := coefficients[j+1];
            end if;

            # We check if the coefficient to be updated is nonzero.
            if function:-GetPrecision(upolys[i]:-upoly[j]) >= m_i+currentJ-j then 
                coefficient += upolys[i]:-upoly[j]:-HomogeneousPart_local(upolys[i]:-upoly[j],
                                                                        m_i+currentJ-j);
            end if;
            # We update the coefficient.
            PowerSeriesObject:-SetHomogeneousPart(upolys[i]:-upoly[j], 
                                                    m_i+currentJ-j, 
                                                    coefficient);
        end do;
    end do;

    # The next iteration to be run.
    state:-currentJ := state:-currentJ + 1;
end proc;

# In charge of knowing when a iteration of MainIteration
# must be run.
local HenselGenerator ::static := proc( self :: PowerSeriesObject,
                                        req_deg :: nonnegint,
                                        asso_up :: UnivariatePolynomialOverPowerSeriesObject,
                                        state :: record,
                                        $)
    local d := DEGREE(state:-init_upop);
    local m_i := self:-GetAncestor(self, ':-m_i');
    local coeff_numb := self:-GetAncestor(self, ':-coeff_numb');

    # Each iteration gives homogeneous components of
    # degree m_i+currentJ.
    # The degree in v, coeff_numb, of each monomial does not change.
    # Only the degree of u is increase, i.e., always goes 
    # from 0 to d-1. 
    # On the other hand, for the currentJ iteration, since the 
    # components are homogeneous of degree m_i+currentJ
    # the degrees of u are gonna be m_i + currentJ - coeff_numb, with 
    # coeff_numb=0,...,d-1.
    while state:-currentJ <= req_deg - (m_i - coeff_numb) do
      asso_up:-MainIteration(state);
    end do;
     
    return self:-HomogeneousPart_local(self, req_deg);

end proc;

# Convert an expanded polynomial p to an array up to degree deg_v.  
local convert_hpoly_from_poly ::static := proc(p::algebraic,
                                                deg_v :: nonnegint,
                                                u::name, 
                                                $)

    ASSERT(p = frontend(expand, [p]));
    if p=0 then
        return Array(0 .. deg_v, 0);
    end if;

    local terms := p; # p is already expanded!
    local ld := degree(terms, u);
    terms := convert(terms, ':-list', ':-`+`');
    if not type(terms, ':-list'(':-polynom')) then
        error "unsupported term: %1", remove(type, terms, ':-polynom')[1];
    else
        terms := ListTools:-Classify(degree, terms, u);
    end if;
    
    return Array(0 .. max(ld, deg_v), map(eq -> lhs(eq) = add(rhs(eq)),  {indices(terms, 'pairs')}));
end proc;

# Extended Hensel construction generator. 
local ehc_gen ::static := proc(self :: PowerSeriesObject, 
                    d :: nonnegint,
                    $)
    return HenselGenerator(self, d, PowerSeriesObject:-GetAncestor(self, ':-asso_up'),
    						PowerSeriesObject:-GetAncestor(self, ':-state'));
end proc;

# Extended Hensel construction generator (multivariate case). 
local ehc_gen_multivariate ::static := proc(self :: PowerSeriesObject, 
                    d :: nonnegint,
                    $)
    local my_pso := PowerSeriesObject:-GetAncestor(self, ':-dummy_pso');
    local state := my_pso:-GetAncestor(my_pso, ':-state');   

    local var;
    local init_vars := state:-init_upop:-Variables(state:-init_upop);
    local v := state:-v;
    local u := state:-u;
    local vars := init_vars minus {v};
    local vars_as_conts := [u=1, seq(var(0)=var, var in vars)];

    return  subs(vars_as_conts, my_pso:-HomogeneousPart(my_pso, d));
end proc;

# Creates a Hensel UPoP from init_fact_in.
local HenselUpop::static := proc(self::UnivariatePolynomialOverPowerSeriesObject, 
									init_fact_in::algebraic,
                                    my_root::algebraic,
                                    state::record, multivariate::truefalse:=false
                                    , $)
    local v := state:-v;
    local u := state:-u;
    local dummy;
    local init_fact := AUTO_EXPAND(dummy, init_fact_in);

	local d := degree(init_fact, v);
    # The output.
    local p_out := Object(UnivariatePolynomialOverPowerSeriesObject, Array(0 .. d), 
    						self:-vname, ':-checkvariables' = false);
    local b := p_out:-upoly;
    local t; 

    # Coefficients as a polynomial in v.
    local coefficients := PolynomialTools:-CoefficientList(init_fact, v);
	for local i from 1 to numelems(coefficients)-1 do 
		local coef := coefficients[i];
        local homo_poly_array := convert_hpoly_from_poly(coef, d, u);
        local deg := numelems(homo_poly_array) - 1;

        # state contains the common information for 
        # all the factors of the Hensel construction.
        # We need to update all the factors at the same
        # time. Hensel, all the Pso coefficients of all
        # the factors.
        b[i-1] := Object(PowerSeriesObject, homo_poly_array, deg,
                           ehc_gen, {u}, 
                           [':-asso_up'=p_out,
                           ':-init_fact'=init_fact,
                           ':-my_root'=my_root, 
                           ':-m_i'=d,
                           ':-coeff_numb'=i-1,
                           ':-state'=state]);
    end do;


    # The factors must be monic.
    local one_gen := proc(self :: PowerSeriesObject, 
                            d :: nonnegint,
                            $)
        return ifelse(d = 0, 1, 0);
    end proc;
    b[numelems(coefficients)-1] :=  Object(PowerSeriesObject, Array(0..0, 1), 0,
                           one_gen, {}, 
                           [':-asso_up'=p_out,
                           ':-init_fact'=init_fact, 
                           ':-my_root'=my_root, 
                           ':-m_i'=d,
                           ':-coeff_numb'=numelems(coefficients)-1,
                           ':-state'=state]);

    return p_out;
end proc;

local multivariate_wrapper::static := proc(self::UnivariatePolynomialOverPowerSeriesObject,
                                            state::record, $)
    local v := state:-v;
    local u := state:-u;
    local var;
    local init_vars := state:-init_upop:-Variables(state:-init_upop);
    local vars := init_vars minus {v};
    local vars_as_conts := [u=1, seq(var(0)=var, var in vars)];
    local coeffs := self:-upoly;
    local new_coeffs := Array(0..DEGREE(self));

    # Multivariate wrapper.
    for local i from 1 to numelems(coeffs)-1 do 
        local c0 := PowerSeriesObject:-HomogeneousPart(coeffs[i-1], 0);
        local c0_subs := subs(vars_as_conts, c0);

        new_coeffs[i-1] := Object(PowerSeriesObject, Array(0..0, c0_subs), 0,
                       ehc_gen_multivariate, vars, 
                       [':-dummy_pso'=coeffs[i-1]]);
    end do;

    new_coeffs[numelems(coeffs)-1] := coeffs[numelems(coeffs)-1];

    return Object(UnivariatePolynomialOverPowerSeriesObject, new_coeffs, 
                            self:-vname, ':-checkvariables' = false);;
end proc;

# To convert the input into a 
# UPoP with PusO coefficients, and then we
# need to apply eval((1/lc)^deg(self,x)*self, v=v*lc) 
# to each self in hensel_polys.
local ConvertToHenselUPoP::static := proc(self_in :: UnivariatePolynomialOverPowerSeriesObject, 
                                            lc::PowerSeriesObject, 
                                            state::record,
                                            mp::list(`=`( name, {:-`*`( { name, name^rational } ), name, name^rational})) := [], 
                                            $)
    local multivariate := state:-multivariate;
    local self;
    if multivariate then 
        self := multivariate_wrapper(self_in, state);
    else 
        self := self_in; 
    end if ;

    local A := map(PuiseuxSeries, self:-upoly, mp);

    local lc_constant_part := lc:-HomogeneousPart(lc, 0);

    # The lc of the input polynomial is different than 1.
    if lc_constant_part <> 1 or lc:-GetAnalyticExpression(lc) <> 1 then
        local lc_as_puso := PuiseuxSeries(lc);
        local lc_inv := PuiseuxSeriesObject:-Inverse(lc_as_puso);

        local d := DEGREE(self);

        for local i from 0 to d-1 do 
            local lc_power := PuiseuxSeriesObject:-Exponentiate(lc_inv, d-i);
            A[i] := PuiseuxSeriesObject:-BinaryMultiply(A[i], lc_power);
        end do;
        A[d] := A[d];
    end if;


    return MultivariatePowerSeries:-UnivariatePolynomialOverPowerSeries(A, self:-vname);
end proc;

# To compute the discriminant of an UPoP up to the degree bnd.
# If mode=analytic, then we first attempt to use the analytic
# expression of self (as long as this is possible). 
# If the above is not possible or if mode=bound, we truncate 
# the coefficients of self up to degree bound, and then compute
# the discriminant. 
export Discrim::static := proc(self::UnivariatePolynomialOverPowerSeriesObject,
                                bnd:: integer := -1, 
                                {mode :: identical(bound, analytic) := ':-bound'},
                                $)
    local self_poly;

    if mode = analytic then 
        self_poly := self:-GetAnalyticExpression(self);
        if type(self_poly, ':-polynom') and self_poly <> undefined then 
            return discrim(self_poly, self:-vname);
        end if;
    end if;

    if bnd < 0 then 
        self_poly := self:-Truncate(self);
    else 
        self_poly := self:-Truncate(self, bnd);
    end if;

    return discrim(self_poly, self:-vname);
end proc;

# All the branches of the polynomial self given by the extended
# Hensel construction.
# self must be monic and square-free.
# For the multivariate case, we need to introduce a new variable,
# say t, then we evaluate all the variables u (except the main 
# variable v) in u=t*u. Now, we can see self_in as a polynomial in the 
# variables t and v, with coefficients in the other variables.
# Then, we can apply the EHC for the bivariate case.
local ExtendedHenselConstructionUnivariate::static := proc(self_in::UnivariatePolynomialOverPowerSeriesObject,
										Hensel_bound :: {integer, identical(FAIL)} := FAIL, 
                                        {returnleadingcoefficient :: {truefalse, identical(automatic)}
                                        := ':-automatic'}, 
                                        {multiplicity::truefalse := false}, $)
	local my_hensel_bound;

	# We check bounds.
    if Hensel_bound <> FAIL then
        my_hensel_bound := Hensel_bound;
    elif self_in:-Hensel_bound <> undefined then
        my_hensel_bound := self_in:-Hensel_bound;
    else
        my_hensel_bound := self_in:-Hensel_bound_static;
    end if;

    # We check that the input polynomial is square-free.
    if self_in:-Discrim(self_in, my_hensel_bound, ':-mode'=analytic)=0 then 
        error "%1 must be square-free", self_in;
    end if;

    # self must be monic.
    local d :=  DEGREE(self_in);
    local lc := self_in:-upoly[d];
    local v := self_in:-MainVariable(self_in);
    local self, i;

    # Note: If lc<>1, then the following change of variables
    # is applied. self := lc^(d-1)*eval(self_in, v=v/lc).
    # Note, that this is equivalent to  process  used in the 
    # for loop. Furthermore, self is now a monic UPoP.
    if lc:-GetAnalyticExpression(lc)=1 then
        self := self_in;
    elif lc:-ApproximatelyZero(lc, my_hensel_bound) then
        error "leading coefficient of %1 may be zero", self_in;
    else 
        local my_upoly := Array(0..d);
        for i from 0 to d-2 do 
            my_upoly[i] := lc^(d-1-i)*self_in:-upoly[i];
        end do;
        my_upoly[d-1] := self_in:-upoly[d-1];
        my_upoly[d] := PowerSeriesObject:-One();
        self := Object(UnivariatePolynomialOverPowerSeriesObject, my_upoly, 
                            v, ':-checkvariables' = false);
    end if;

    local vars := self:-Variables(self) minus {v};

    # Multivariate case.
    # The univariate case is written in such a way that it allows an
    # UPoP with multivariate coefficients as an input. This, however, is 
    # not used in practice . The reason is that it is better to convert first 
    # multivariate UPoP to the case univariate polynomial with univariate 
    # coefficients and then called this algorithm.
	local u, var; 
    local multivariate:= false;
    if numelems(vars)>1 then 
        u := t;
        multivariate := true;
    else 
        u := op(vars);
    end if;

    # Monic initial factors of the construction.
    local d_hat, delta_hat, Facs;
    d_hat, delta_hat, Facs := self:-initialFactors(self, my_hensel_bound, v, u, ':-multivariate'=multivariate);
	local size := numelems(Facs);

	# We lift our univariate factors.
    local factors:= [seq(Groebner:-Homogenize(Facs[i], u), i=1..size)];
	local my_roots:= [seq(-1*coeff(Facs[i][1], v, 0), i=1..size)];
    local init_factors := [seq(eval(factors[i][1]^factors[i][2], u=u^(delta_hat)), i=1..size)];

    # We compute the Lagrange interpolation polynomials.
    # Note: Facs are univariate polynomials.
    local lag_polys;
    lag_polys := self:-FastLagrangeInterpolationPolynomials(Facs, v, u, delta_hat);
 
	# The output.
    # Note: if lc<>1 then we need to revert the change of variable 
    # applied at the beginning. Thus, self_in must be equal to
    # (1/lc)^(d-1)*prod(eval(f, v=v*lc), f in hensel_polys). This is equivalent
    # to lc*(1/lc)^d*prod(eval(f, v=v*lc), f in hensel_polys), which can be 
    # rewritten as:
    # lc*prod(eval((1/lc)^deg(f,x)*f, v=v*lc), f in hensel_polys).
    local hensel_polys; 
    local lc_included := true;
    if returnleadingcoefficient<>false then
        local lc_constant_part := lc:-HomogeneousPart(lc, 0);
        if lc_constant_part = 1 and lc:-GetAnalyticExpression(lc) = 1 then
            if returnleadingcoefficient=true then 
                hensel_polys := Array(0..size);
                hensel_polys[0] := Object(UnivariatePolynomialOverPowerSeriesObject, 
                                            Array(0..0, [PuiseuxSeries(1)]), v);
                
            else 
                hensel_polys := Array(1..size);
                lc_included := false;
            end if;
        else
            local lc_upop := Object(UnivariatePolynomialOverPowerSeriesObject, [PuiseuxSeries(lc)], 
                                            v, ':-checkvariables' = false);
            hensel_polys := Array(0..size);
            hensel_polys[0] := lc_upop;
        end if;
    else 
        hensel_polys := Array(1..size);
        lc_included := false;
    end if;

    # We initialize the fast_product and the LiftedFactors for the
    # FastProduct algorithm.
    # Case k=0.
    local fast_product := Array(1..size+1);
    local LiftedFactors := Array(1..size);
    for local j from 1 to size do
        fast_product[j] :=   Array(0..-1, 0);
        LiftedFactors[j] := Array(0..0, [init_factors[j]]);
    end do;

    # Trick: We pass hensel_factors as a name at initialization,
    # and after that we hensel_factors := hensel_polys.
    local hensel_factors;
    # The common record of all the hensel_polys.
    local state := Record('currentJ'=1, 'upolys'=hensel_factors, ':-lag_polys'=lag_polys,
                            'init_upop'=self, ':-d_hat'=d_hat, ':-delta_hat'=delta_hat,
                            ':-lc_included'=lc_included, ':-fast_product'=fast_product,
                            ':-LiftedFactors'=LiftedFactors, ':-multivariate'=multivariate,
                            ':-multiplicity'=multiplicity, ':-v'=v, ':-u'=u);

    # We create new UPoPs objects.
	for i from 1 to size do 
    	hensel_polys[i] := HenselUpop(self, init_factors[i], my_roots[i], state, multivariate);
    end do;
    hensel_factors := hensel_polys;

    local poly, mp;
    if multivariate then 
        mp := [seq(var=var, var in vars)];
    else 
        mp := [u=u^(1/d_hat)];
    end if;
    # First, we need to convert each factor into a 
    # UPoP with PusO coefficients, and then we
    # need to apply eval((1/lc)^deg(poly,x)*poly, v=v*lc) 
    # to each poly in hensel_polys.
    local my_output;
    if lc_included then 
	   my_output := [[hensel_polys[0]], seq([poly, lc, state, mp], poly in hensel_polys[1..])];
    else 
       my_output := [seq([poly, lc, state, mp], poly in hensel_polys[1..])];
    end if;

    return my_output;
end proc;

## Here the code of the classic version is finished.
## We now have the code for multiplicity of the newton polynomial
## and the "multivariate" case.

# Extended Hensel construction generator (multivariate case). 
# To see a power series in x1,..,xn as a power series in a
# variable t such that the change of variables xi=xi*t is 
# applied.
local substitute_gen ::static := proc(self :: PowerSeriesObject, 
                    d :: nonnegint,
                    $)
    local dummy_pso := PowerSeriesObject:-GetAncestor(self, ':-dummy_pso');
    local u := PowerSeriesObject:-GetAncestor(self, ':-u');
    local eqns := PowerSeriesObject:-GetAncestor(self, ':-eqns');
    local init_vars := dummy_pso:-Variables(dummy_pso);

    return  subs(eqns, dummy_pso:-HomogeneousPart(dummy_pso, d));
end proc;

# To see a upop with coefficients in x1,..,xn as a upop in a with
# coefficients in the variable t such that the change of variables xi=xi*t is 
# applied to each coefficient.
local multivariate_substitute :: static := proc(self :: UnivariatePolynomialOverPowerSeriesObject,
                                                u::name,
                                                eqns :: {set, list, thistype}(name = {PowerSeriesObject, algebraic}),
                                                 $)
    
    # We first try to use the analytic expression of self.
    local ana := self:-GetAnalyticExpression(self);
    local v := self:-MainVariable(self);
    if type(ana, ':-polynom') and ana <> undefined then 
        local new_ana := subs(eqns, ana);

        return UnivariatePolynomialOverPowerSeries(new_ana, v);
    end if;

    # Non polynomial analytic expression.
    local p;
    local my_upoly := Array(0..DEGREE(self));
    local old_upoly := self:-upoly;

    for local i from 0 to DEGREE(self) do
        local c0 := PowerSeriesObject:-HomogeneousPart(old_upoly[i], 0);
        local c0_subs := subs(eqns, c0);
        my_upoly[i] := Object(PowerSeriesObject, Array(0..0, c0_subs), 0,
                            substitute_gen, {u}, 
                            [':-dummy_pso'=old_upoly[i], ':-u'=u, ':-eqns'=eqns]);
    end do;

    return Object(UnivariatePolynomialOverPowerSeriesObject, my_upoly, 
                            v, ':-checkvariables' = false);
end proc;

# When the newton polynomial contains factors with multiplicity, this algorithm
# is called in a recursive way to apply the ExtendedHenselConstructionUnivariate
# to each of this factors.
local ExtendedHenselConstructionMultiplicity::static := proc(self_in :: UnivariatePolynomialOverPowerSeriesObject, 
                                            lc::PowerSeriesObject, 
                                            state::record,
                                            mp::list(`=`( name, {:-`*`( { name, name^rational } ), name, name^rational})) := [], 
                                            Hensel_bound :: {integer, identical(FAIL)} := FAIL, 
                                            $)
    # We apply a TschirnhausenTransformation to change the root of self
    # (in other case, self would not be square free).
    local delta_hat := state:-delta_hat;
    local d_hat := state:-d_hat;
    local n := DEGREE(self_in);
    local init_fact := PowerSeriesObject:-GetAncestor(self_in:-upoly[n], ':-init_fact');
    local u := state:-u;
    local my_root := PowerSeriesObject:-GetAncestor(self_in:-upoly[n], ':-my_root');
    local a := PowerSeries(-1*my_root*u^delta_hat);

    # We apply the transformation x+my_root*u^delta_hat.
    local self := self_in:-TschirnhausenTransformation(self_in, a, 1, powerseriesmode=true);

    # self must be univariate already, so we use ':-returnleadingcoefficient'=false.
    local my_factors := self:-ExtendedHenselConstructionUnivariate(self, Hensel_bound,
                                                     ':-returnleadingcoefficient'=false,
                                                     ':-multiplicity'=true);

    # mp must be updated with the info of the previously computed d_hat in each factor.
    # We also need to apply the inverse transformation.
    for local i from 1 to numelems(my_factors) do 
        local eq;
        my_factors[i][1] := self:-TschirnhausenTransformation(my_factors[i][1], a, 1, detransform=true, powerseriesmode=true);
        my_factors[i][4] := [seq(lhs(eq)=rhs(eq)^(1/d_hat), eq in my_factors[i][4])];
    end do;
    
    # Recursive call.
    local f;
    my_factors := [seq(ifelse(DEGREE(f[1])>1, 
                            self:-ExtendedHenselConstructionMultiplicity(op(f), Hensel_bound), f), f in my_factors)];


    return op(my_factors);
end proc;

# General algorithm to wrap all the cases.
export ExtendedHenselConstruction::static := proc(self_in::UnivariatePolynomialOverPowerSeriesObject,
                                        Hensel_bound :: {integer, identical(FAIL)} := FAIL, 
                                        t::name:= NULL,
                                        {returnleadingcoefficient :: {truefalse, identical(automatic)}
                                       := ':-automatic'}, 
                                       {output :: {thistype,list}(identical(factorization, changeofvariables)) := factorization},
                                        $)

    local v := self_in:-MainVariable(self_in);
    local vars := self_in:-Variables(self_in) minus {v};
    
    local u, eqns, var; 
    local multivariate:= false;
    local self;

    if numelems(vars)>1 then 
        # Multivariate case: We apply the change of variables
        # ui=ui(0)*t where ui is vars. Then, we apply the 
        # ExtendedHenselConstructionUnivariate normally. 
        # Note: here we see now the variables ui as constants.
        if t = NULL then 
            error "in multivariate mode, ExtendedHenselConstruction requires a variable name";
        end;

        u := t;
        multivariate := true;
        eqns := [seq(var=var(0)*t, var in vars)];
        self := self_in:-multivariate_substitute(self_in, u, eqns);
    else 
        # Univariate case.
        u := op(vars);
        self := self_in;
    end if;

    # Univariate call.
    local my_factors := self:-ExtendedHenselConstructionUnivariate(self, Hensel_bound, _options['returnleadingcoefficient']);

    # Recursive call for the multiplicity case.
    local f;
    local new_factors := [seq(ifelse(DEGREE(f[1])>1, 
                            self:-ExtendedHenselConstructionMultiplicity(op(f), Hensel_bound), f), f in my_factors)];

    # We convert the output to UPoPs with Puiseux series coefficients.
    if multivariate and output=':-changeofvariables' then 
        eqns := [t=1, seq(var(0)=var, var in vars)];
        return [[seq(ifelse(numelems(f)>1, self:-ConvertToHenselUPoP(op(f)), op(f)), f in new_factors)], eqns];
    else    
        return [seq(ifelse(numelems(f)>1, self:-ConvertToHenselUPoP(op(f)), op(f)), f in new_factors)];
    end if;
end proc;

# To set a Hensel_bound or a Hensel_bound_static.
export SetHenselBound ::static := proc(_self :: UnivariatePolynomialOverPowerSeriesObject,
                                                 bound::nonnegint,
                                                 {instance :: truefalse := true},
                                                 $)
    local old_value;
    
    if instance then
        old_value := _self:-Hensel_bound;
        _self:-Hensel_bound := bound;
    else
        old_value := _self:-Hensel_bound_static;
        _self:-Hensel_bound_static := bound;
    end if;

    return old_value;
end proc;