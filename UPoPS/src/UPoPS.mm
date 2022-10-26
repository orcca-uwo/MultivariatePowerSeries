
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
    dstyle :: UPOPSDS_TYPE,
    # Bound for the Puiseux theorem
    Puiseux_bound_static :: static(nonnegint) := 10,
    Puiseux_bound :: {nonnegint, identical(undefined)} := undefined; 

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
        local mystyle, displayresult;

        if user_dstyle <> [] then
            mystyle := user_dstyle;
        elif self:-dstyle <> [] then
            mystyle := self:-dstyle;
        else
            mystyle := defaultDisplayStyle;
        end if;
        local psstyle := select(type, mystyle, 'PSDS_ENTRYTYPE');
        
        local min_deg := min(DEGREE(self), eval(':-maxdegree', [op(mystyle), ':-maxdegree'=infinity]));
        local nterms := eval(':-maxterms', [op(mystyle), ':-maxterms'=50]);
        local result := Array(1 .. 0);
        # _Z shouldn't be documented
        local x := ifelse(self:-vname <> undefined, self:-vname, '_Z');

        displayresult := function:-Display(self:-upoly[0], psstyle, output, 
                                           ':-updateterms' = 'nterms');

        if output = "string" then
            result ,= sprintf("(%s)", displayresult);
        else
            result ,= T:-mfenced(T:-mrow(displayresult));
        end if;

        local d;
        for d from 1 to min_deg while nterms > 0 do
            displayresult := function:-Display(
                self:-upoly[d], [':-maxterms' = nterms, op(psstyle)], output, ':-updateterms' = 'nterms');

            if output = "string" then
                result ,= sprintf(" + (%s) %a", displayresult, x^d);
            else
                result ,= T:-mn("&plus;"), T:-mfenced(T:-mrow(displayresult)), T:-mo("&InvisibleTimes;"), apply(T:-Typeset, x^d);
            end if;
        end do;

        # The last degree included, at least partially, is d-1.
        if d <= DEGREE(self) then
            if output = "string" then
                result ,= " + ...";
            else
                result ,= T:-mn("&plus;"), T:-mn("&hellip;");
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

$include "MultivariatePowerSeries/UPoPS/src/basic_routines.mm"
$include "MultivariatePowerSeries/UPoPS/src/basic_arithmetic.mm"
$include "MultivariatePowerSeries/UPoPS/src/weierstrass_preparation.mm"
$include "MultivariatePowerSeries/UPoPS/src/factorization.mm"
$include "MultivariatePowerSeries/UPoPS/src/Puiseux_factorization.mm"

$undef AUTO_EXPAND
$undef DEGREE

end module:
