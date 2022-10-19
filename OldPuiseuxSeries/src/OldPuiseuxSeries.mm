local module OldPuiseuxSeriesObject()
option object;

local 
	pso, # Power Series object
	expDivs::table, #Table mapping each variable to the
							 # Denominator of its exponent 
    dstyle::PSDS_TYPE; #a display style for this object


local ModuleApply::static := proc()
	return Object(OldPuiseuxSeriesObject, _passed);
end proc;


	#constructor of the class
local ModuleCopy::static := proc(new :: OldPuiseuxSeriesObject, 
								 old :: OldPuiseuxSeriesObject,
								 pso::PowerSeriesObject,
								 expDivs::{table, list(name=posint)} := [],
                                 dstyle::PSDS_TYPE := [],
                                 $)
	local l;
	new:-pso := pso;
    new:-dstyle := dstyle;

	if type(expDivs, ':-table') then
		l := [entries(expDivs, ':-nolist', ':-pairs')];
		if not type(l, ':-list'(':-name'=':-posint')) then
			error("invalid input: %1 expects its %-2 argument, %3, to consist of entries of type %4, but received %5", 
									procname, 4, ':-expDivs', ':-name' =':-posint', expDivs);
		end if;
	else
  		l := expDivs;
	end if; 

	l := select((eqn, vars) -> lhs(eqn) in vars, l, pso:-Variables(pso));
	new:-expDivs := table(':-sparse' = 1, l);
end proc;		

# SetDisplayStyle
# to set the dstyle of the input object
export SetDisplayStyle :: static:= proc(_self :: OldPuiseuxSeriesObject, s :: PSDS_TYPE, $)
    _self:-dstyle := s; 
end proc;

export Display :: static := proc(_self :: OldPuiseuxSeriesObject,
                                 user_dstyle :: list := [],
                                 output :: identical("typeset", "string") := "typeset",
                                 $)
uses T = Typesetting;
local mystyle;
    if user_dstyle <> [] then
        mystyle := user_dstyle;
    elif _self:-dstyle <> [] then
        mystyle := _self:-dstyle;
    else
        mystyle := PowerSeriesObject:-GetDefaultDisplayStyle();
    end if;

local pair;
    return _self:-pso:-Display(
        _self:-pso, mystyle, ifelse(output = "typeset", "full", "fullstring"),
        ':-objectname' = "OldPuiseuxSeries", ':-substitutions' = {
            seq(lhs(pair) = lhs(pair)^(1/rhs(pair)), pair in [indices(_self:-expDivs, ':-pairs')])});
end proc;


# ModulePrint: print _self :: OldPuiseuxSeriesObject 
local ModulePrint :: static := proc(_self :: OldPuiseuxSeriesObject, $)
    if IsWorksheetInterface() then
        return _self:-Display(_self);
    else
        return convert(_self:-Display(_self, [], "string"), ':-name');
    end if;
end proc;

# Module Deconstruct
local ModuleDeconstruct :: static := proc(_self :: OldPuiseuxSeriesObject, $)
    return 'Object'('OldPuiseuxSeriesObject', _self:-pso, _self:-expDivs);
end proc;	 

# Basic routines
$include "MultivariatePowerSeries/OldPuiseuxSeries/src/getsAndSets.mm"
$include "MultivariatePowerSeries/OldPuiseuxSeries/src/makeCompatible.mm"
# Arithmetic 
$include "MultivariatePowerSeries/OldPuiseuxSeries/src/addition.mm"
$include "MultivariatePowerSeries/OldPuiseuxSeries/src/multiplication.mm"
$include "MultivariatePowerSeries/OldPuiseuxSeries/src/equality.mm"

end module; 
