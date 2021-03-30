unprotect( 'MultivariatePowerSeries' ):

MultivariatePowerSeries := module()
    option package;

    local   
        PowerSeriesObject, # PSObject
        UnivariatePolynomialOverPowerSeriesObject; #UPoPSObject

$define PSDS_ENTRYTYPE \
    :-identical(:-maxterms, :-precision) = {:-nonnegint, :-identical(:-infinity )}

$define PSDS_TYPE \
    list(PSDS_ENTRYTYPE)

$define UPOPSDS_ENTRYTYPE \
    {PSDS_ENTRYTYPE, :-identical(:-maxdegree) = {:-nonnegint, :-identical(:-infinity)}}

$define UPOPSDS_TYPE \
    list(UPOPSDS_ENTRYTYPE)

    # Set the "global" default display style in PSObject and UPoPSObject 
    export
        SetDefaultDisplayStyle := proc(s :: UPOPSDS_TYPE, $)
            PowerSeriesObject:-SetDefaultDisplayStyle(select(type, s, 'PSDS_ENTRYTYPE')); 
            UnivariatePolynomialOverPowerSeriesObject:-SetDefaultDisplayStyle(s);
        end proc;

    # Set the "local" attribute, dstyle, in PSObject and UPoPSObject 
    export    
        SetDisplayStyle := proc(obj :: {PowerSeriesObject, UnivariatePolynomialOverPowerSeriesObject}, 
                                s :: {PSDS_TYPE, UPOPSDS_TYPE}, 
                                $)
            obj:-SetDisplayStyle(obj, s);
        end proc;

    # a wrapper for PSO:-DeepCopy and UPoPS:-DeepCopy
    export 
        Copy := proc(obj :: {PowerSeriesObject, UnivariatePolynomialOverPowerSeriesObject}, $)
            return obj:-DeepCopy(_passed);
        end proc;

    # Display method for PSO and UPoPS
    export Display ::static := proc(obj :: {PowerSeriesObject, UnivariatePolynomialOverPowerSeriesObject}, 
                                    displayStyle :: UPOPSDS_TYPE := [],
                                    $)
            obj:-Display(obj, displayStyle);
        end proc;

    # a wrapper for PSO:-HomogenousPart 
    export 
        HomogeneousPart := proc(obj :: PowerSeriesObject, 
                                d :: nonnegint, 
                                $)
                return obj:-HomogeneousPart(_passed);
        end proc;

    # a wrapper for PSO:-GetCoefficient and UPoPS:-GetCoefficient
    export 
        GetCoefficient := proc(obj :: {PowerSeriesObject, UnivariatePolynomialOverPowerSeriesObject}, 
                            idx :: {nonnegint, polynom},
                            $)
            obj:-GetCoefficient(_passed);
        end proc;

    # a wrapper for PSO:-Precision
    export 
        Precision := proc(obj :: PowerSeriesObject, $)
            return obj:-GetPrecision(obj);
        end proc;

    # a wrapper for UPoPS:-Degree
    export 
        Degree := proc(obj :: UnivariatePolynomialOverPowerSeriesObject, $)
            return obj:-Degree(_passed);
        end proc;

    # a wrapper for UPoPS:-MainVariable
    export 
        MainVariable := proc(obj :: UnivariatePolynomialOverPowerSeriesObject, $)
            return obj:-MainVariable(_passed);
        end proc;

    # a wrapper for PS:-GetAnalyticExpression
    export 
        GetAnalyticExpression := proc(obj :: {PowerSeriesObject, UnivariatePolynomialOverPowerSeriesObject}, $)
            return obj:-GetAnalyticExpression(_passed);
        end proc;

    # a wrapper for UPoPS:-UpdatePrecision
    export 
        UpdatePrecision := proc(obj :: {PowerSeriesObject, UnivariatePolynomialOverPowerSeriesObject}, 
                                prec :: nonnegint,
                                $)
            return obj:-UpdatePrecision(_passed);
        end proc;

    # PowerSeries
    # like a constructor for PowerSeriesObject
    # in fact, a wrapper for PSO:-One, PSO:-Zero, PSO:-Constant, PSO:-FromPolynomial, PSO:-FromProcedure, and UPoPS:-ConvertToPowerSeries
    # TODO: support rational function? 
    export
        PowerSeries := proc(p :: {UnivariatePolynomialOverPowerSeriesObject, PowerSeriesObject, procedure, appliable_module, polynom, numeric, algnum, abstract_rootof, radical}, 
                            {variables::set(name) := undefined},
                            {analytic::algebraic := undefined},
                            {check::truefalse := true},
                            {expand::truefalse := true},
                            $)
            if type(p, {procedure, appliable_module}) and 
            not type(p, {UnivariatePolynomialOverPowerSeriesObject, PowerSeriesObject}) then
                return PowerSeriesObject:-FromProcedure(p, analytic, variables, _options['check', 'expand']);
            end if;
            if analytic <> undefined or check <> true or expand <> true or variables <> undefined then 
                error "invalid input: the options variables, analytic, expand, and check cannot be "
                "used with this calling sequence";
            end if;
            if type(p, UnivariatePolynomialOverPowerSeriesObject) then 
                return UnivariatePolynomialOverPowerSeriesObject:-ConvertToPowerSeries(p);
            elif type(p, PowerSeriesObject) then
                return p:-DeepCopy(p);
            elif type(p, numeric) then
                return PowerSeriesObject:-Constant(p);
            elif type(p, {numeric, algnum, algnum^fraction}) then
                local expanded := Algebraic:-Expand(p);
                return PowerSeriesObject:-Constant(expanded);
            else 
                return PowerSeriesObject:-FromPolynomial(p);
            end if;
        end proc;

    # UnivariatePolynomialOverPowerSeries
    # like a constructor for UnivariatePolynomialOverPowerSeriesObject
    # in fact, a wrapper for UPoPS:-Zero, UPoPS:-One, UPoPS:-Constant, UPoPS:-FromPolynomial, UPoPS:-FromPowerSeriesList
    # TODO: support rational function? 
    export 
        UnivariatePolynomialOverPowerSeries := proc(p :: {polynom, algnum, algnum^fraction,
                                                          {list, Array, Vector, thistype}(PowerSeriesObject),
                                                          UnivariatePolynomialOverPowerSeriesObject},
                                                    x :: name,
                                                    $)
            if type(p, {':-numeric', ':-algnum', ':-algnum'^':-fraction'}) then
                local pp := ifelse(type(p, ':-numeric'), p, Algebraic:-Expand(p));
                return UnivariatePolynomialOverPowerSeriesObject:-Constant(pp);
            elif type(p, ':-polynom') then
                return UnivariatePolynomialOverPowerSeriesObject:-FromPolynomial(p, x);
            elif type(p, PowerSeriesObject) then
                return p:-ConvertToUPoPS(p, x);
            elif type(p, UnivariatePolynomialOverPowerSeriesObject) then
                if _npassed = 1 or p:-MainVariable(p) in {undefined, x} then
                    return p:-DeepCopy(p);
                else
                    error "you specified %1 as the main variable, but the main variable of the first "
                    "argument is %2", x, p:-MainVariable(p);
                end if;
            else
                local pp := (if type(p, {':-list', ':-And'(':-Array', [0] &under [lowerbound])}) then
                                 p;
                             elif rtable_num_dims(p) <> 1 then
                                 error "invalid input: expected 1-dimensional rtable, but found dimension %1",
                                 rtable_num_dims(p);
                             else
                                 ArrayTools:-Alias(p, [0 .. upperbound(p) - lowerbound(p)]);
                             end if);
                
                if ormap(pso -> x in pso:-Variables(pso), pp) then
                    error "expected list or Array of power series objects independent of %1, but found "
                    "some that were dependent on it", x;
                end if;
                
                return Object(UnivariatePolynomialOverPowerSeriesObject, pp, x);
            end if;
        end proc;

    # a wrapper for PSO:-Truncate and UPoPS:-Truncate
    export 
        Truncate := proc(obj :: {PowerSeriesObject, UnivariatePolynomialOverPowerSeriesObject})
            return obj:-Truncate(_passed);
        end proc;

    # a wrapper for PSO:-ApproximatelyEqual and UPoPS:-ApproximatelyEqual
    export 
        ApproximatelyEqual := proc(obj1 :: {PowerSeriesObject, UnivariatePolynomialOverPowerSeriesObject}, 
                        obj2 :: depends(obj1),
                        deg :: integer := -2,
                        $)
                return obj1:-ApproximatelyEqual(_passed);
        end proc;

    # a wrapper for PSO:-ApproximatelyZero and UPoPS:-ApproximatelyZero
    export 
        ApproximatelyZero := proc(obj :: {PowerSeriesObject, UnivariatePolynomialOverPowerSeriesObject}, 
                    deg :: integer := -2,
                    $)
            obj:-ApproximatelyZero(_passed);
        end proc;

    # a wrapper for PSO:-IsUnit 
    export 
        IsUnit := proc(obj :: {PowerSeriesObject, UnivariatePolynomialOverPowerSeriesObject}, $)
            obj:-IsUnit(_passed);
        end proc;

    # create a geometric series 
    export 
        GeometricSeries := proc(vars :: {name, nonemptylist(name)}, $)
            return PowerSeriesObject:-GeometricSeries(vars);
        end proc;

    # Create a sum_of_all_monomials
    export
        SumOfAllMonomials := proc(vars :: nonemptylist(name), $)
            return PowerSeriesObject:-SumOfAllMonomials(vars);
        end proc;

    # Add
    export 
        Add := proc(obj :: seq({PowerSeriesObject, UnivariatePolynomialOverPowerSeriesObject, polynom, algebraic}),
                    {coefficients :: list({polynom, algebraic}) := []},
                    $)
            local objs := select(type, [obj], ':-object'); 
            if numelems(objs) = 0 then 
                error "invalid input: expected a sequence of inputs containing objects, but no object provided";
            elif membertype(UnivariatePolynomialOverPowerSeriesObject, objs) then
                UnivariatePolynomialOverPowerSeriesObject:-NaryAdd(_passed);
            elif membertype(PowerSeriesObject, objs) then
                PowerSeriesObject:-NaryAdd(_passed);
            else 
                error "invalid input: expected a sequence of inputs containing power series objects "
                "and/or univariate polynomial over power series objects, but found %1", objs;
            end if;
        end proc;

    # Subtract
    export 
        Subtract := proc(obj1 :: {PowerSeriesObject, UnivariatePolynomialOverPowerSeriesObject},
                        obj2 :: depends(obj1),
                        $)
            return obj1:-BinarySub(_passed);
        end proc;

    # Negate 
    export 
        Negate := proc(obj :: {PowerSeriesObject, UnivariatePolynomialOverPowerSeriesObject}, $)
            obj:-Negate(_passed);
        end proc;

    # Multiply
    export 
        Multiply := proc(obj :: seq({PowerSeriesObject, UnivariatePolynomialOverPowerSeriesObject, polynom, algebraic}), $)
            local objs := select(type, [obj], ':-object'); 
            if numelems(objs) = 0 then 
                error "invalid input: expected a sequence of inputs containing objects, but no object provided";
            elif type(objs, 'list'(PowerSeriesObject)) then
                PowerSeriesObject:-NaryMultiply(_passed);
            elif type(objs, 'list'(UnivariatePolynomialOverPowerSeriesObject)) then
                UnivariatePolynomialOverPowerSeriesObject:-NaryMultiply(_passed);
            else 
                error "invalid input: expected a list of PowerSeriesObjects or a list of \
                        UnivariatePolynomialOverPowerSeriesObjects, but received %1", objs;
            end if;
        end proc;

    # Exact Quotient PowerSeries
    export 
        Divide := proc(obj1 :: PowerSeriesObject,
                            obj2 :: PowerSeriesObject,
                            $)
            return obj1:-BinaryExactQuotient(_passed);
        end proc;

    # Inverse PowerSeries
    export 
        Inverse := proc(obj :: PowerSeriesObject, $)
            obj:-Inverse(_passed);
        end proc;

    # Exponentiate 
    export 
        Exponentiate := proc(obj :: {PowerSeriesObject, UnivariatePolynomialOverPowerSeriesObject}, 
                        n :: integer,
                        $)
            obj:-Exponentiate(_passed);
        end proc;
    
    # TaylorShift UPoPS
    export 
        TaylorShift := proc(obj :: UnivariatePolynomialOverPowerSeriesObject, 
                            c :: {numeric, algebraic, algnum, algnumext, abstract_rootof},
                            $)
            obj:-TaylorShift(_passed);
        end proc;

    # Weierstrass Preparation UPoPS
    export 
        WeierstrassPreparation := proc(obj :: UnivariatePolynomialOverPowerSeriesObject, $)
            return obj:-WeierstrassPreparation(obj);
        end proc;

    # EvaluateAtOrigin UPoPS
    export 
        EvaluateAtOrigin := proc(obj :: UnivariatePolynomialOverPowerSeriesObject, $)
            return obj:-EvaluateAtOrigin(_passed);
        end proc;

    # HenselFactorize UPoPS
    export 
        HenselFactorize := proc(obj :: UnivariatePolynomialOverPowerSeriesObject,
                                {returnleadingcoefficient :: {truefalse, identical(automatic)}
                                 := ':-automatic'},
                                $)
            return obj:-HenselFactorize(_passed);
        end proc;

    export
        Variables := proc(obj :: {UnivariatePolynomialOverPowerSeriesObject, PowerSeriesObject}, $)
            return obj:-Variables(obj);
        end proc;

$include "PowerSeries/src/PowerSeries.mm"
$include "UPoPS/src/UPoPS.mm"

$undef PSDS_ENTRYTYPE 
$undef PSDS_TYPE 
$undef UPOPSDS_ENTRYTYPE
$undef UPOPSDS_TYPE

end module:

protect( 'MultivariatePowerSeries' ):
#savelib( 'MultivariatePowerSeries' ):
 
