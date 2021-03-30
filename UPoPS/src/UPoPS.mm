
module UnivariatePolynomialOverPowerSeriesObject()
option object;

local 
    # upoly : A Maple Array to store PowerSeriesObjects and represent a dense univariate polynomial 
    #         we don't use Maple list or other data types as we frequently update, copy, 
    #         and in some cases reshape upoly. 
    upoly :: Array,
    # vname : tha name of main variable 
    vname :: name,
    # dstyle : a display style for a specific UPoPS
    dstyle :: UPOPSDS_TYPE; 

# return Expanded of x wehere d is a dummy variable for the sake of efficiency 
$define AUTO_EXPAND(d, x) \
    ifelse(type((d := x), ':-polynom'(':-complex'(':-numeric'))), expand(d), Algebraic:-Expand(subsindets(d, ':-radical', convert, RootOf)))

# return degree of x, i.e. upperbound(x:-upoly) 
$define DEGREE(x) \
    upperbound(x:-upoly)

# ModuleApply : to make UPoPS an object factory
local
    ModuleApply ::static := proc()
        Object(UnivariatePolynomialOverPowerSeriesObject, _passed);
    end proc;

# ModuleCopy: called when a new object of a class is created via a call to Object
local 
    ModuleCopy ::static := proc(new :: UnivariatePolynomialOverPowerSeriesObject,
                                old :: UnivariatePolynomialOverPowerSeriesObject,
                                upoly :: {Array, list},
                                vname :: name := undefined,
                                dstyle :: UPOPSDS_TYPE := [],
                                {checkvariables :: truefalse := true},
                                $)
        if _npassed > 2 then
            if type(upoly, ':-Array') then
                if lowerbound(upoly) <> 0 then
                    error "invalid input: expected Array with lower bound 0, but found %1", lowerbound(upoly);
                end if;
                new:-upoly := upoly;
            else
                new:-upoly := Array(0 .. numelems(upoly)-1, upoly);
            end if;
            new:-dstyle := dstyle;
            new:-vname := vname;
            local pso;
            if checkvariables and vname <> undefined
            and vname in `union`(seq(function:-Variables(pso), pso in new:-upoly)) then
                local i := min(select(j -> vname in function:-Variables(new:-upoly[j]),
                                      [seq(0 .. upperbound(new:-upoly))]));
                error "invalid input: the main variable, %1, occurs in %2", vname, (
                    if i = 0 then
                        "the constant coefficient";
                    elif i = 1 then
                        sprintf("the coefficient of %a", vname);
                    else
                        sprintf("the coefficient of %a^%d", vname, i);
                    end if);
            end if;
        elif type(old:-upoly, 'Array') then
            new:-upoly := copy(old:-upoly, 'deep');
            new:-vname := old:-vname;
            new:-dstyle := old:-dstyle;
        else
            error "you cannot copy the original UnivariatePolynomialOverPowerSeriesObject object";
        end if;
    end proc;

local 
defaultDisplayStyle ::static := [];

# SetDefaultDisplayStyle
# To set the local global variable defaultDisplayStyle 
export 
    SetDefaultDisplayStyle ::static := proc(s :: UPOPSDS_TYPE, $)
        defaultDisplayStyle := s;
    end proc;

# SetDisplayStyle
# to set the dstyle of the input object
export    
    SetDisplayStyle ::static:= proc(obj :: UnivariatePolynomialOverPowerSeriesObject, s :: UPOPSDS_TYPE, $)
        obj:-dstyle := s; 
    end proc;

# Display
# to disply the input object, 
# user_dstyle : will overwrite the default and self:-dstyle style
# output : (doesn't need to be documented) used for internal purposes
export 
    Display ::static := proc(self :: UnivariatePolynomialOverPowerSeriesObject, 
                             user_dstyle :: list := [],
                             output :: identical("typeset", "string") := "typeset",
                             $)
        uses T = Typesetting;
        local mystyle;
        if user_dstyle <> [] then
            mystyle := user_dstyle;
        elif self:-dstyle <> [] then
            mystyle := self:-dstyle;
        else
            mystyle := defaultDisplayStyle;
        end if;
        
        local min_deg := min(DEGREE(self), eval(':-maxdegree', [op(mystyle), ':-maxdegree'=infinity]));
        local result := Array(0 .. min_deg);
        # _Z shouldn't be documented
        local x := ifelse(self:-vname <> undefined, self:-vname, '_Z');
        local displayresult := PowerSeriesObject:-Display(self:-upoly[0], select(type, mystyle, 'PSDS_ENTRYTYPE'), output);
        if output = "string" then
            result[0] := sprintf("(%s) + ", displayresult);
        else
            result[0] := T:-mfenced(T:-mrow(displayresult)),  T:-mn("&plus;");
        end if;
        # it will update defaultDisplayStyle... 
        for local d from 1 to min_deg - 1 do
            displayresult := PowerSeriesObject:-Display(self:-upoly[d], select(type, mystyle, 'PSDS_ENTRYTYPE'), output);
            if output = "string" then
                result[d] := sprintf("(%s) %a + ", displayresult, x^d);
            else
                result[d] := T:-mfenced(T:-mrow(displayresult)), T:-mo("&InvisibleTimes;"), apply(T:-Typeset, x^d), T:-mn("&plus;");
            end if;
        end do;

        displayresult := PowerSeriesObject:-Display(self:-upoly[min_deg], select(type, mystyle, 'PSDS_ENTRYTYPE'), output);
        if min_deg > 0 then
            if output = "string" then
                result[min_deg] := sprintf("(%s) %a", displayresult, x^min_deg);
            else
                result[min_deg] := T:-mfenced(T:-mrow(displayresult)), T:-mo("&InvisibleTimes;"), apply(T:-Typeset, x^min_deg);
            end if;
        else
            if output = "string" then
                result[min_deg] := sprintf("(%s)", displayresult);
            else
                result[min_deg] := T:-mfenced(T:-mrow(displayresult));
            end if;
        end if;

        if min_deg < DEGREE(self) then
            if output = "string" then
                result[min_deg] := cat(result[min_deg], " + ...");
            else
                result[min_deg] := result[min_deg], T:-mn("&plus;"), T:-mn("&hellip;");
            end if;
        end if;

        if output = "string" then
            return cat("< UnivariatePolynomialOverPowerSeries: ", seq(result), " >");
        else
            return T:-mfenced(T:-mrow(T:-mn("UnivariatePolynomialOverPowerSeries:   "),
                                      seq(result)), ':-open' = "&lsqb;", ':-close' = "&rsqb;");
        end if;
    end proc; 

# ModulePrint: print self :: UnivariatePolynomialOverPowerSeriesObject 
local
    ModulePrint :: static := proc(self :: UnivariatePolynomialOverPowerSeriesObject, $)
        if IsWorksheetInterface() then
            return Display(self);
        else
            return convert(Display(self, [], "string"), ':-name');
        end if;
    end proc;

# Module Deconstruct
local
  ModuleDeconstruct :: static := proc(self :: UnivariatePolynomialOverPowerSeriesObject, $)
    return 'Object'('UnivariatePolynomialOverPowerSeriesObject', self:-upoly);
  end proc;

$include "UPoPS/src/basic_routines.mm"
$include "UPoPS/src/basic_arithmetic.mm"
$include "UPoPS/src/weierstrass_preparation.mm"
$include "UPoPS/src/factorization.mm"

$undef AUTO_EXPAND
$undef DEGREE

end module:
