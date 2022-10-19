# Make compatible two PuSO.
#
# _self have an associate cone given by _sefl:-rays 
# and the same happens with other. This cones represent
# the change of variables applied to each PuSO. If we 
# want to add or multiply two PuSO, they must depend of
# the same change of variables, i.e., they must have the 
# same associate cone. The reason for this is 
# that we want to add or multiply their internal pso.
#
# Our first step is to make our set of rays (the union
# of _self:-rays, other:-rays and extraRays) integer.
# The reason for this is that we base our algorithm 
# in the paper:
# Formal Laurent Series in Several Variables.
# by:
# Ainhoa Aparicio Monforte, Manuel Kauers.
# This document explains how to code Laurent series.
# Thus, we compute the lcm of all the denominator of
# all our rays, and we call it "denominator". Then we 
# multiply all our rays by "denominator". Note that 
# this is equivalent to have all the rays with the 
# same denominator, but we do not write this denominator
# until the end of our computations.
#
# Now, we look for a minimal set of rays, which generated
# cone contains the union of the cone associated to _self
# and the cone associated to other. Then, we write the rays
# of the originals cones in terms of the new cone, i.e., we 
# compute the original changes of variables in terms of 
# the new one. By the properties of our cone, we know these 
# two computations are going to output a change of variables
# with positive integers exponents. Hence, we are able 
# to apply these changes of variables to _self:-pso and
# other:-pso without problems.
#
# A key observation is that all this computations are done
# with integer cones. However, we want to get back our rational
# change of variable. Thus, we multiply by the inverse of 
# "denominator". 
local BinaryMakeCompatible::static := proc(_self :: PuiseuxSeriesObject, 
											other :: PuiseuxSeriesObject,  
											extraRays ::list(list(rational)) :=[], $)
	local newRays, pso1, pso2, cv1, cv2, M, newOrdCV, i, newOrd, newLso,
		  s, r;

	local rays := [op(_self:-rays), op(other:-rays), op(extraRays)];
	local denominator := lcm(seq(seq(denom(s), s in r) , r in rays));

	#########
    # Puiseux Series adjustment.
    if denominator<>1 then
		rays := denominator*rays;
	end if;
	#########

	newOrd := _self:-ord;
		
	# We get a new set of rays that generate a cone C
	# such that _self and other can be seen as elements of K_C[[X]]
	newRays := _self:-MakeRaysCompatible(_self, rays, newOrd);

	# We get a new order for the internal power series objects
	newOrdCV := _self:-ComputeNewOrdCV(_self, numelems(newRays), 
									   _self:-ordCV, other:-ordCV);
	
	M := Matrix(newRays);

	# We make a change of variables in the internal pso using the newRays
	# and the newOrdCV
	pso1 := _self:-MakePsoCompatible(_self, M, 
									   newOrdCV, denominator);
	pso2 := other:-MakePsoCompatible(other, M, 
									   newOrdCV, denominator);

	newRays := map(convert, newRays, ':-list');

	#########
    # Puiseux Series adjustment.
    if denominator<>1 then
	newRays := (1/denominator)*newRays;
	end if;
	#########

	return newOrd, newOrdCV, newRays, pso1, pso2;
end proc;

local ComputeNewOrdCV::static := proc(_self :: PuiseuxSeriesObject,
									  num_newRays::integer,
									  ordCV1::list(name),
									  ordCV2::list(name) := [],
									  $)
# We get a new order for the internal power series objects
	local i;
	local newOrdCV := [op(ordCV1), op(ordCV2)];

	newOrdCV := ListTools:-MakeUnique(newOrdCV);

	if numelems(newOrdCV) > num_newRays then
		newOrdCV := newOrdCV[1 .. num_newRays];
	elif numelems(newOrdCV) < num_newRays then
		newOrdCV := [op(newOrdCV), seq(cat(':-_W', i), 
									i = 1 .. num_newRays - numelems(newOrdCV))];
	end if;

	return newOrdCV;
end proc;

# We look for a minimal set of rays that generates a cone, which contains
# the cone generated by the rays in "rays".
# The grevlex order can be defined by a sequence of weight vectors. So,
# we use this fact to compare our rays and compute a minimal set of
# rays called e, such that the cone generated by e, contains the cone
# generated by "rays". 
local MakeRaysCompatible::static := proc(_self :: PuiseuxSeriesObject, 
										 rays::list(list(rational)), 
										 ord::list(name),  $) 
	local i, II, v, m, functionals, scale, vscaled;

	local d1 := nops(ord);
	local d2 := nops(rays);

	local e := Array(1 .. 0); 

	local grevlexFunctional := _self:-GrevlexFunctional;

	# We create a matrix using the rays
	local w := Matrix(rays);

	for m from 1 to d1 do
		# We select all the columns of w such that the sum
		# of its first m-coordinates is greater than 0
		II := select(i -> grevlexFunctional(m, w[i]) > 0, [seq(1..d2)]);

		# We "normalize" the columns selected in II
		if nops(II) = 1 then
			v := [seq(w[i] , i in II)];
		elif II = [] then
			next;
		else 
			functionals := [seq(grevlexFunctional(m, w[i]), i in II)];
			scale := igcd(functionals);
			v := [seq(w[i], i in II)];
			v := [seq(v[i] / (functionals[i]/scale), i=1..numelems(v))];
		end if;

		# We look for the grevlex smallest column of v
		v := _self:-Smallest(_self,[seq(convert(v[i], ':-list'), i=1..nops(II))]);
		# We look for a integer vector that's grevlex-less than ray
		v:= _self:-LookForGreatestGrevlexLess(_self, v);

		# We save v
		e ,= Vector[':-row'](v);
		vscaled := e[:-upperbound(e)] / GrevlexFunctional(m, v);
		for i in II do
			w[i] := w[i] -GrevlexFunctional(m, w[i]) * vscaled;
		end do;
		
	end do;

	return [seq(convert(v, ':-list'), v in e)]; 
end proc;

# We make an appropriate change of variables in pso to see _self as
# an element of the cone generate by M
local MakePsoCompatible::static := proc(_self :: PuiseuxSeriesObject,
										M::Matrix,  
										ordCV::list(name), 
										denominator::integer,
										$) 
	
	local r;
	local rays := map(convert, _self:-rays, ':-Vector');
	#########
    # Puiseux Series adjustment.
    if denominator<>1 then
		rays := denominator*rays;
	end if;

	# We do Linear Algebra to compute the change of variables in terms
	# of the new variables and rays.
	local result := [seq(LinearAlgebra:-LinearSolve(M^+, r), r in rays)];

	local cv := [seq(mul(ordCV^~r), r in result)];
	cv := _self:-ordCV=~cv;

	return _self:-pso:-Substitute(cv, _self:-pso);
end proc;

# We make the orders from _self and other compatible, if possible
local MakeOrdCompatible::static := proc(_self :: PuiseuxSeriesObject,
										ord1::list(name),
										ord2::list(name),   
										 $) 
	
	local S1, S2, j1, j2, v;
	local new_ord := Array(1..0);
	
	# We compute the intersection of the orders as sets
	local inter := convert(ord1, ':-set') intersect convert(ord2, ':-set');

	# We compute the complement of inter
	local A := convert(ord1, ':-set') union convert(ord2, ':-set') minus inter;
	
	# We remove from ord1 and ord2 the elements of A.
	# If the lists generated are different, then it means the orders are
	# not compatible
	if remove(`in`,ord1, A) <> remove(`in`,ord2, A) then
		error "the order of Puiseux series %1 and %2 are not compatible", ord1, ord2;
	end if;

	# We generate the new order
	j1 := 1;
	j2 := 1;
	for v in inter do
		S1, j1 := _self:-SublistLessThanV(_self, ord1, v, j1);
		S2, j2 := _self:-SublistLessThanV(_self, ord2, v, j2);

		new_ord ,= _self:-DisjointListUnion(_self, S1, S2), v;
	end do;

	new_ord ,= _self:-DisjointListUnion(_self, ord1[j1..], ord2[j2..]);

	return convert(new_ord, ':-list');	
end proc;

local SublistLessThanV := proc(_self :: PuiseuxSeriesObject,
										ls :: list(name),
										v::name,
										i::integer,   
										 $) 

	local j := ListTools:-Search(v, ls);
	return [op(i..j-1,ls)], j+1;

end proc;

local DisjointListUnion  :: static := proc(_self :: PuiseuxSeriesObject,
    s1 :: list,
    s2 :: list,
    $)
    local i1 := 1;
    local i2 := 1;
    local result := Array(1..0);

    while i1 <= numelems(s1) and i2 <= numelems(s2) do
        if {s1[i1], s2[i2]}[1] = s1[i1] then
            result ,= s1[i1++];
        else
            result ,= s2[i2++];
        end if;
    end do;

    if i1 <= numelems(s1) then
        result ,= op(i1.., s1);
    elif i2 <= numelems(s2) then
        result ,= op(i2.., s2);
    end if;

    return seq(result);
end proc;

# We extend a lso with respect to a compatible order.
# We take the rays in _self and we inject them in the 
# bigger space generated by new_ord.
local ExtendPuiseuxSeriesObject::static := proc(_self :: PuiseuxSeriesObject,
										        new_ord::list(name),   
										 $) 
	local p, i, r, v;

	local ord_as_set := convert(_self:-ord, ':-set');
	local ray := [seq(ifelse(v in ord_as_set, v, 0), v in new_ord)];
	local new_rays := [seq(eval(ray, _self:-ord=~r), r in _self:-rays )];

	return Object(PuiseuxSeriesObject, _self:-pso, new_ord, _self:-ordCV, new_rays, _self:-e);
end proc;

# To convert a PuSO to PSO whenever is possible.
export ConvertToPSO ::static := proc(_self :: PuiseuxSeriesObject, $)
    # We get the analytic expression on _self.
    local p :=  _self:-GetAnalyticExpression(_self);
    p := normal(p);

    # We check first if p is a number, a polynomial or a rational
    # polynomial.
    if type(p, 'numeric') then
        return PowerSeriesObject:-Constant(p);
    elif type(p, ':-polynom') and p<>undefined then
        return PowerSeriesObject:-FromPolynomial(p);
    elif type(p, 'ratpoly') and p<>undefined then
        local num := PowerSeriesObject:-FromPolynomial(numer(p));
        local den := PowerSeriesObject:-FromPolynomial(denom(p));
        return num / den;
    else
    	# Finally, we try to create a PSO using a procedure. 
	    local vars := convert(_self:-ord, ':-set');
	    local ps := Object(PowerSeriesObject, Array(0 .. 0), 0, PowerSeriesObject:-ConvertToPSO_gen,
	    			 vars, ["A" = _self, "cache_table" = Array(0..-1), "is_error"=false],
	    			 ':-initializeconstantterm');
    	
    	return ps; 
    end if;

end proc;

# We apply a change of variables in a univariate PuSO.
# This is an auxiliary function of the Puiseux theorem for UPoPS.
export ChangeOfVariablesForPuiseuxTheorem::static := proc(_self::PuiseuxSeriesObject, 
                                      q::nonnegint, p::integer, k::integer, $)
    local var:= _self:-ord;

    # We check that we have a univariate PuSO.
    # If we have a constant PuSO, there is nothing
    # to do.
    if numelems(var) > 1 then 
		error "invalid input: %1 must be a"
              " polynomial in one variable", _self;
    elif numelems(var) = 0 then
    	return _self:-ConvertToPSO(_self);
	end if;   

	local e:= copy(_self:-e);
	var := var[1];
	e[var] := q*_self:-e[var];

	local rays := [];
	if _self:-rays<>[] then
		rays := [[op(_self:-rays[1])*q]];
	end if;
                              
	local puso := Object(PuiseuxSeriesObject, _self:-pso, _self:-ord, 
									 		  _self:-ordCV, rays, e);

	puso := puso*(1/var)^(k*p);

	return puso:-ConvertToPSO(puso);
end proc; 