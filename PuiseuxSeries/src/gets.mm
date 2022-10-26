#Gets.

# Get pso.
export GetPowerSeries::static := proc(_self :: PuiseuxSeriesObject, $)
	return _self:-pso; 
end proc;

# Get order of the Puiseux series variable.
export GetPuiseuxSeriesOrder::static := proc(_self :: PuiseuxSeriesObject, $)
	return _self:-ord; 
end proc;

# Get order of the power series variable.
export GetPowerSeriesOrder::static := proc(_self :: PuiseuxSeriesObject, $)
	return _self:-ordCV; 
end proc;

# Get the monomial that multiplies the PuSO.
export GetMonomial::static := proc(_self :: PuiseuxSeriesObject, $)
	local x;

	return mul(x^_self:-e[x], x in _self:-ord); 
end proc;

# Get the rays.
export GetRays::static := proc(_self :: PuiseuxSeriesObject, $)
	return _self:-rays; 
end proc;

# Get the analytic expression.
export GetAnalyticExpression::static := proc(_self :: PuiseuxSeriesObject, $)
	local ana := _self:-pso:-GetAnalyticExpression(_self:-pso);

	if ana=undefined then
		return undefined;
	else
		local mon := _self:-GetMonomial(_self);
		local change_of_variables := _self:-ChangeOfVariables(_self);

		return mon*subs(change_of_variables, ana);
	end if;
end proc;

# Get the order of a Puiseux series. For a Puiseux series
# in the form P=sum_{m=M}^{inf} a_m X^{m/n}, for some M \in Z 
# (not necessarily >= 0). The the order of P is defined as
# min{m/n : a_m <> 0}.
export GetOrder::static := proc(_self :: PuiseuxSeriesObject, 
								bnd :: nonnegint := FAIL, $)
	local my_Puiseux_bound;

	if numelems(_self:-ord) > 1 then 
		error "invalid input: %1 must be a"
              " polynomial in one variable", _self;
	end if;

	# We check the Puiseux's bound.
    if bnd <> FAIL then
        my_Puiseux_bound := bnd;
    elif _self:-Puiseux_bound <> undefined then
        my_Puiseux_bound := _self:-Puiseux_bound;
    else
        my_Puiseux_bound := _self:-Puiseux_bound_static;
    end if;

	local ana := _self:-GetAnalyticExpression(_self);
	if ana=0 then
		return infinity;
	end if;

	if type(ana,numeric) then
		return 0;
	end if;

	if _self:-ord = [] then
		return 0;
	end if;

	local p;

	for local i from 0 to my_Puiseux_bound do
		p := function:-HomogeneousPart(_self:-pso, i);
		p := subs(_self:-ChangeOfVariables(_self), p);

		if p<>0 then
			if _self:-rays <> [] then
				return i*_self:-rays[1][1] + _self:-e[_self:-ord[1]];
			else 
				return _self:-e[_self:-ord[1]];
			end if;
		end if;
	end do;

	error "the order of %1 cannot be determined"
          " with the current bound %2", _self, my_Puiseux_bound;
end proc;