#Gets

# Get pso
export GetPowerSeries::static := proc(_self :: LaurentSeriesObject, $)
	return _self:-pso; 
end proc;

# Get order of the Laurent series variable
export GetLaurentSeriesOrder::static := proc(_self :: LaurentSeriesObject, $)
	return _self:-ord; 
end proc;

# Get order of the power series variable
export GetPowerSeriesOrder::static := proc(_self :: LaurentSeriesObject, $)
	return _self:-ordCV; 
end proc;

# Get the monomial that multiplies the LSO
export GetMonomial::static := proc(_self :: LaurentSeriesObject, $)
	local x;

	return mul(x^_self:-e[x], x in _self:-ord); 
end proc;

# Get the rays
export GetRays::static := proc(_self :: LaurentSeriesObject, $)
	return _self:-rays; 
end proc;

# Get the rays
export GetAnalyticExpression::static := proc(_self :: LaurentSeriesObject, $)
	local ana := _self:-pso:-GetAnalyticExpression(_self:-pso);

	if ana=undefined then
		return undefined;
	else
		local mon := _self:-GetMonomial(_self);
		local change_of_variables := _self:-ChangeOfVariables(_self);

		return mon*subs(change_of_variables, ana);
	end if;
end proc;