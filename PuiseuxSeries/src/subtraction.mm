# BinarySub
# to subtract two objects
export BinarySub ::static := proc(_self :: PuiseuxSeriesObject,
                                  other :: PuiseuxSeriesObject, $)
   local neg_other := other:-Negate(other);

   return _self:-BinaryAdd(_self, neg_other);
end proc;

# Negate
# to negate the input object
export Negate ::static := proc(_self :: PuiseuxSeriesObject, $)
  local pso := _self:-pso:-Negate(_self:-pso);
  
  return Object(PuiseuxSeriesObject, pso, _self:-ord, _self:-ordCV, 
                                          _self:-rays, _self:-e);
end proc;
