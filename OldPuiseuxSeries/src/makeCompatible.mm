# To make to Puiseux object compatible. That means: getting the same denominator in the exponents
# of the variables that are in _self and other
# We return PowerSeriesObjects
local BinaryMakeCompatible :: static 
              := proc(_self :: {OldPuiseuxSeriesObject, PowerSeriesObject
                              , UnivariatePolynomialOverPowerSeriesObject},
                      other :: {OldPuiseuxSeriesObject, PowerSeriesObject
                              , UnivariatePolynomialOverPowerSeriesObject},
                      $)
  local new_expDivs, A, B, C;

  if type(_self, OldPuiseuxSeriesObject) then
    if type(other, PowerSeriesObject) then
      return BinaryMakeCompatibleWithPSO(_self, other);
    elif type(other, OldPuiseuxSeriesObject) then
      return BinaryMakeCompatibleWithPuSO(_self, other);
    else 
      C := other:-ConvertToPowerSeries(other);

      return BinaryMakeCompatibleWithPSO(_self, C);
    end if;
  elif type(_self, PowerSeriesObject) then
    if type(other, OldPuiseuxSeriesObject) then
      new_expDivs, B, A := BinaryMakeCompatibleWithPSO(other, _self);

      return new_expDivs, A, B;
    elif type(other, PowerSeriesObject) then
      return table(':-sparse' = 1, []), _self, other;
    else 
      C := other:-ConvertToPowerSeries(other);

      return table(':-sparse' = 1, []), _self, C;
    end if;
  elif type(_self, UnivariatePolynomialOverPowerSeriesObject) then
    if type(other, OldPuiseuxSeriesObject) then
      C := _self:-ConvertToPowerSeries(_self);
      new_expDivs, B, A := BinaryMakeCompatibleWithPSO(other, C);

      return new_expDivs, A, B;
    elif type(other, PowerSeriesObject) then
      C := _self:-ConvertToPowerSeries(_self);

      return table(':-sparse' = 1, []), C, other;
    else 
      A := _self:-ConvertToPowerSeries(_self);
      B := other:-ConvertToPowerSeries(other);

      return table(':-sparse' = 1, []), A, B;
    end if;
  end if;
                                                 
end proc;
 

local BinaryMakeCompatibleWithPSO :: static := proc(_self :: OldPuiseuxSeriesObject,
            other :: PowerSeriesObject,
            $)
  local v, v_B, B;

  # We make change of variables to look for equivalent fractions in the exponents with the same
  # denominator
  v_B := {seq(v = v^(_self:-expDivs[v]), v in other:-Variables(other))};
  B := other:-Substitute(v_B, other);

  return _self:-expDivs, _self:-pso, B;
                                        
end proc;

local BinaryMakeCompatibleWithPuSO :: static := proc(_self :: OldPuiseuxSeriesObject,
            other :: OldPuiseuxSeriesObject,
            $)
  local v, new_expDivs, v_A, v_B, A, B;

  # We get the lcm between the denominators of each exponent
  new_expDivs := table([seq(v=ilcm(_self:-expDivs[v], other:-expDivs[v]), 
                                                    v in (_self:-Variables(_self) 
                                                    union other:-Variables(other)))]);

  # We make change of variables to look for equivalent fractions in the exponents with the same
  # denominator
  v_A := {seq(v = v^(new_expDivs[v]/_self:-expDivs[v]), v in _self:-Variables(_self))};
  A := _self:-pso:-Substitute(v_A, _self:-pso);

  v_B := {seq(v = v^(new_expDivs[v]/other:-expDivs[v]), v in other:-Variables(other))};
  B := other:-pso:-Substitute(v_B, other:-pso);

  return new_expDivs, A, B;
                                             
end proc;

