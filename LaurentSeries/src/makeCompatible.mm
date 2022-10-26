# Make compatible two LSO
local BinaryMakeCompatible::static := proc(_self :: LaurentSeriesObject, 
											other :: LaurentSeriesObject,  
											extraRays ::list(list(integer)) :=[], $)
	local newRays, pso1, pso2, cv1, cv2, M, newOrdCV, i, newOrd, newLso;

	local rays := [op(_self:-rays), op(other:-rays), op(extraRays)];

	local newSelf := _self;
	local newOther := other;

	# We check if the Laurent Series orders of _self and other agrees
	if _self:-ord <> other:-ord then
  		# If the ord of _self and other is compatible, we extend 
		# other and _self to the Laurent Series space generated 
		# by this new order.
		newOrd := _self:-MakeOrdCompatible(_self, _self:-ord, other:-ord);
		
		newOther := other:-ExtendLaurentSeriesObject(other, newOrd);
		newSelf := _self:-ExtendLaurentSeriesObject(_self, newOrd);
	end if;

	newOrd := newSelf:-ord;
		
	# If one of then is a constant, we handle it separetely
	if newSelf:-ordCV = [] or newOther:-ordCV = [] then
		newRays := [op(newSelf:-rays), op(newOther:-rays)];
		newOrdCV := [op(newSelf:-ordCV), op(newOther:-ordCV)];
		pso1 := newSelf:-pso;
		pso2 := newOther:-pso;
	else
		# We get a new set of rays that generate a cone C
		# such that newSelf and newOther can be seen as elements of K_C[[X]]
		newRays := newSelf:-MakeRaysCompatible(newSelf, rays, newOrd);

		# We get a new order for the internal power series objects
		newOrdCV := _self:-ComputeNewOrdCV(_self, numelems(newRays), 
										   newSelf:-ordCV, newOther:-ordCV);
		
		M := Matrix(newRays);

		# We make a change of variables in the internal pso using the newRays
		# and the newOrdCV
		pso1 := newSelf:-MakePsoCompatible(newSelf, M, newOrdCV);
		pso2 := newOther:-MakePsoCompatible(newOther, M, newOrdCV);

		newRays := map(convert, newRays, ':-list');
	end if;

	return newOrd, newOrdCV, newRays, pso1, pso2;
end proc;

local ComputeNewOrdCV::static := proc(_self :: LaurentSeriesObject,
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

# We look for a minimal set of rays that generates a cones which contains
# the cone generate by the rays of _self and other
local MakeRaysCompatible::static := proc(_self :: LaurentSeriesObject, 
										 rays::list(list(integer)), 
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
			v := [seq(v[i] / (functionals[i] / scale), i=1..numelems(v))];
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

# We make an appropiate change of variables in pso to see _self as
# an element of the cone generate by M
local MakePsoCompatible::static := proc(_self :: LaurentSeriesObject,
										M::Matrix,  
										ordCV::list(name), $) 
	
	local r;
	local rays := map(convert, _self:-rays, ':-Vector');

	# We do Linear Algebra to compute the change of variables in terms
	# of the new variables and rays.
	local result := [seq(LinearAlgebra:-LinearSolve(M^+, r), r in rays)];

	local cv := [seq(mul(ordCV^~r), r in result)];
	cv := _self:-ordCV=~cv;

	return _self:-pso:-Substitute(cv, _self:-pso);
end proc;

# We make the orders from _self and other compatible, if possible
local MakeOrdCompatible::static := proc(_self :: LaurentSeriesObject,
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
		error "the order of Laurent series %1 and %2 are not compatible", ord1, ord2;
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

local SublistLessThanV := proc(_self :: LaurentSeriesObject,
										ls :: list(name),
										v::name,
										i::integer,   
										 $) 

	local j := ListTools:-Search(v, ls);
	return [op(i..j-1,ls)], j+1;

end proc;

local DisjointListUnion  :: static := proc(_self :: LaurentSeriesObject,
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
# bigger space generated by new_ord
local ExtendLaurentSeriesObject::static := proc(_self :: LaurentSeriesObject,
										        new_ord::list(name),   
										 $) 
	local p, i, r, v;

	local ord_as_set := convert(_self:-ord, ':-set');
	local ray := [seq(ifelse(v in ord_as_set, v, 0), v in new_ord)];
	local new_rays := [seq(eval(ray, _self:-ord=~r), r in _self:-rays )];

	return Object(LaurentSeriesObject, _self:-pso, new_ord, _self:-ordCV, new_rays, _self:-e);
end proc;
