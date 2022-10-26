# Get pso
export GetPso::static := proc(_self :: OldPuiseuxSeriesObject, $)
		return _self:-pso; 
end proc;

# Get pso
export GetExpDivs::static := proc(_self :: OldPuiseuxSeriesObject, $)
		return eval(_self:-expDivs); 
end proc;

# Get Variables
export Variables::static := proc(_self :: OldPuiseuxSeriesObject, $)
		return _self:-pso:-Variables(_self:-pso); 
end proc;
