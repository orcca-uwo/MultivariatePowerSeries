unprotect( 'MultivariatePowerSeries' ):

MultivariatePowerSeries := module()
    option package;

# The order of the types is important: checking type complexcons is potentially SLOOOWWWWW
$define COEFFICIENT_TYPE \
    Or(numeric, radnumext, algnumext, complexcons)
    
    local   
        PowerSeriesObject, # PSObject
        UnivariatePolynomialOverPowerSeriesObject, #UPoPSObject
        PuiseuxSeriesObject, #Puiseux series object
        OdeManipulator; #Sub-package OdeManipulator 

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
        SetDisplayStyle := proc(obj :: {PowerSeriesObject, 
                                        PuiseuxSeriesObject,
                                        UnivariatePolynomialOverPowerSeriesObject}, 
                                s :: {PSDS_TYPE, UPOPSDS_TYPE}, 
                                $)
            obj:-SetDisplayStyle(obj, s);
        end proc;

    # a wrapper for PSO:-DeepCopy and UPoPS:-DeepCopy
    export 
        Copy := proc(obj :: {PowerSeriesObject, 
                             UnivariatePolynomialOverPowerSeriesObject,
                             PuiseuxSeriesObject}, $)
            return obj:-DeepCopy(_passed);
        end proc;

    # Display method for PSO and UPoPS
    export Display ::static := proc(obj :: {PowerSeriesObject, 
                                            PuiseuxSeriesObject,
                                            UnivariatePolynomialOverPowerSeriesObject}, 
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
        GetAnalyticExpression := proc(obj :: {PowerSeriesObject, 
                                              UnivariatePolynomialOverPowerSeriesObject,
                                              PuiseuxSeriesObject}, $)
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
    export
        PowerSeries := proc(p :: {UnivariatePolynomialOverPowerSeriesObject, PowerSeriesObject,
                                  procedure, appliable_module, polynom, COEFFICIENT_TYPE, ratpoly, algebraic}, 
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
            elif type(p, 'numeric') then
                return PowerSeriesObject:-Constant(p);
            elif type(p, 'COEFFICIENT_TYPE') then
                local expanded := Algebraic:-Expand(p);
                return PowerSeriesObject:-Constant(expanded);
            elif type(p, 'polynom') then
                return PowerSeriesObject:-FromPolynomial(p);
            elif type(p, 'ratpoly') then
                local num := PowerSeriesObject:-FromPolynomial(numer(p));
                local den := PowerSeriesObject:-FromPolynomial(denom(p));
                return num / den;
            else
                local method_option := ifelse(assigned(_EnvPSMethod), ':-method' = _EnvPSMethod, NULL);
                return PowerSeriesObject:-FromAlgebraicExpression(p, method_option);
            end if;
        end proc;

    # UnivariatePolynomialOverPowerSeries
    # like a constructor for UnivariatePolynomialOverPowerSeriesObject
    # in fact, a wrapper for UPoPS:-Zero, UPoPS:-One, UPoPS:-Constant, UPoPS:-FromPolynomial, UPoPS:-FromPowerSeriesList
    export 
        UnivariatePolynomialOverPowerSeries := proc(p :: {polynom, COEFFICIENT_TYPE, ratpoly, PowerSeriesObject,
                                                         {list, Array, Vector}({PowerSeriesObject, PuiseuxSeriesObject}),
                                                          UnivariatePolynomialOverPowerSeriesObject},
                                                    x :: name,
                                                    $)
            local pp;

            if type(p, 'COEFFICIENT_TYPE') then
                pp := ifelse(type(p, ':-numeric'), p, Algebraic:-Expand(p));
                return UnivariatePolynomialOverPowerSeriesObject:-Constant(pp);
            elif type(p, ':-polynom') then
                return UnivariatePolynomialOverPowerSeriesObject:-FromPolynomial(p, x);
            elif type(p, ':-ratpoly') then
                local den := denom(p);
                local den_inv;
                if has(den, x) then
                    error "invalid input: you specified %1 as the main variable, but it occurs in the "
                    "denominator, %2", x, den;
                end if;

                try
                    den_inv := PowerSeriesObject:-FromPolynomial(den)^(-1);

                    return thisproc(numer(p), x) * den_inv;
                catch "not invertible":
                    local p_as_puso := Object(PuiseuxSeriesObject, p);

                    return PuiseuxSeriesObject:-ConvertToUPoPS(p_as_puso, x);
                    
                end try;

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
                if ormap(type, p, PuiseuxSeriesObject) then
                    pp := map(PuiseuxSeries, p);
                else 
                    pp := p;
                end if;

                pp := (if type(pp, {':-list', ':-And'(':-Array', [0] &under [lowerbound])}) then
                                 pp;
                     elif rtable_num_dims(pp) <> 1 then
                         error "invalid input: expected 1-dimensional rtable, but found dimension %1",
                         rtable_num_dims(pp);
                     else
                         ArrayTools:-Alias(pp, [0 .. upperbound(p) - lowerbound(p)]);
                     end if);
                
                if ormap(pso -> x in pso:-Variables(pso), pp) then
                    error "expected list or Array of power/Puiseux series objects independent of %1, but found "
                    "some that were dependent on it", x;
                end if;
                
                return Object(UnivariatePolynomialOverPowerSeriesObject, pp, x);
            end if;
        end proc;

    export 
        UnivariatePolynomialOverPuiseuxSeries := UnivariatePolynomialOverPowerSeries;

    # a wrapper for PSO:-Truncate and UPoPS:-Truncate
    export 
        Truncate := proc(obj :: {PowerSeriesObject, 
                                 UnivariatePolynomialOverPowerSeriesObject,
                                 PuiseuxSeriesObject})
            return obj:-Truncate(_passed);
        end proc;

    # a wrapper for PSO:-ApproximatelyEqual and UPoPS:-ApproximatelyEqual
    export 
        ApproximatelyEqual := proc(obj1 :: {PowerSeriesObject, 
                                            UnivariatePolynomialOverPowerSeriesObject,
                                            PuiseuxSeriesObject}, 
                        obj2 :: depends(obj1),
                        deg :: integer := -2,
                        {force :: truefalse := false},
                        {mode :: identical(:-powerseries, :-absolute) := ':-powerseries'},
                        $)
                return obj1:-ApproximatelyEqual(_passed);
        end proc;

    # a wrapper for PSO:-ApproximatelyZero and UPoPS:-ApproximatelyZero
    export 
        ApproximatelyZero := proc(obj :: {PowerSeriesObject, 
                                          UnivariatePolynomialOverPowerSeriesObject,
                                          PuiseuxSeriesObject}, 
                    deg :: integer := -2,
                    {force :: truefalse := false},
                    {mode :: identical(:-powerseries, :-absolute) := ':-powerseries'},
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
        Add := proc(obj :: seq({PowerSeriesObject, 
                                UnivariatePolynomialOverPowerSeriesObject, 
                                PuiseuxSeriesObject,
                                polynom, algebraic}),
                    {coefficients :: list({polynom, algebraic}) := []},
                    $)
            local objs := select(type, [obj], ':-object'); 
            
            if numelems(objs) = 0 then 
                error "invalid input: expected a sequence of inputs containing objects, but no object provided";
            elif membertype(UnivariatePolynomialOverPowerSeriesObject, objs) then
                UnivariatePolynomialOverPowerSeriesObject:-NaryAdd(_passed);
            elif membertype(PuiseuxSeriesObject, objs) then
                if coefficients = [] then
                    PuiseuxSeriesObject:-NaryAdd(obj);
                else 
                    error "invalid input: expected a only sequence of inputs containing Puiseux series objects "
                ", but found %1", coefficients;
                end if;
            elif membertype(PowerSeriesObject, objs) then
                PowerSeriesObject:-NaryAdd(_passed);
            else 
                error "invalid input: expected a sequence of inputs containing power series objects "
                "and/or univariate polynomial over power series objects and/or Puiseux series objects,"
                " but found %1", objs;
            end if;
        end proc;

    # Subtract
    export 
        Subtract := proc(obj1 :: {PowerSeriesObject, PuiseuxSeriesObject,
                                  UnivariatePolynomialOverPowerSeriesObject},
                        obj2 :: depends(obj1),
                        $)
            return obj1:-BinarySub(_passed);
        end proc;

    # Negate 
    export 
        Negate := proc(obj :: {PowerSeriesObject, PuiseuxSeriesObject,
                               UnivariatePolynomialOverPowerSeriesObject}, $)
            obj:-Negate(_passed);
        end proc;

    # Multiply
    export 
        Multiply := proc(obj :: seq({PowerSeriesObject, 
                                     UnivariatePolynomialOverPowerSeriesObject, 
                                     PuiseuxSeriesObject,
                                     polynom, algebraic}), $)
            local objs := select(type, [obj], ':-object'); 

            if numelems(objs) = 0 then 
                error "invalid input: expected a sequence of inputs containing objects, but no object provided";
            elif membertype(UnivariatePolynomialOverPowerSeriesObject, objs) then
                UnivariatePolynomialOverPowerSeriesObject:-NaryMultiply(_passed);
            elif membertype(PuiseuxSeriesObject, objs) then
                PuiseuxSeriesObject:-NaryMultiply(_passed);
            elif membertype(PowerSeriesObject, objs) then
                PowerSeriesObject:-NaryMultiply(_passed);
            else 
                error "invalid input: expected a sequence of inputs containing power series objects "
                "and/or univariate polynomial over power series objects and/or Puiseux series objects,"
                " but found %1", objs;
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
        Inverse := proc(obj :: {PowerSeriesObject, PuiseuxSeriesObject}, 
                        pso_bound :: nonnegint := FAIL,
                        extra_bound :: nonnegint := FAIL,
                        {useanalytic :: boolean := true},$)
            obj:-Inverse(_passed);

        end proc;

    # Substitute PowerSeries
    export
        Substitute := proc(eqns, obj :: PowerSeriesObject, $)
            obj:-Substitute(_passed);
        end proc;
    
    # Exponentiate 
    export 
        Exponentiate := proc(obj :: {PowerSeriesObject, PuiseuxSeriesObject,
                        UnivariatePolynomialOverPowerSeriesObject}, 
                        n :: integer,
                        $)
            obj:-Exponentiate(_passed);
        end proc;
    
    # TaylorShift
    export 
        TaylorShift := proc(obj :: {PowerSeriesObject, UnivariatePolynomialOverPowerSeriesObject})
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

    # Extended Hensel contruction UPoPS.
    export ExtendedHenselConstruction := proc(obj :: UnivariatePolynomialOverPowerSeriesObject)
        return obj:-ExtendedHenselConstruction(_passed); 
    end proc;

    # To set a Hensel_bound or a Hensel_bound_static.
    export SetHenselBound := proc(bound::nonnegint,
                                   obj :: UnivariatePolynomialOverPowerSeriesObject := NULL,
                                   $)
        if obj=NULL then
            return function:-SetHenselBound(UnivariatePolynomialOverPowerSeriesObject,
                                        bound, 'instance' = false);
        else
            return function:-SetHenselBound(obj, bound); 
        end if;
    end proc;

    export
        Variables := proc(obj :: {UnivariatePolynomialOverPowerSeriesObject, 
                                  PowerSeriesObject, PuiseuxSeriesObject}, $)
            return obj:-Variables(obj);
        end proc;

    # Tschirnhausen Transformation of an univariate polynomial with 
    # Puiseux series coefficients.
    export TschirnhausenTransformation := proc(obj :: UnivariatePolynomialOverPowerSeriesObject)
        return obj:-TschirnhausenTransformation(_passed); 
    end proc;

    # Puiseux Facorization of an univariate polynomial with 
    # Puiseux series coefficients.
    export PuiseuxFactorize := proc(obj :: UnivariatePolynomialOverPowerSeriesObject)
        return obj:-PuiseuxTheorem(_passed); 
    end proc;

    # Puiseux series exports.
    export PuiseuxSeries := proc(puso)
        if type(puso, PuiseuxSeriesObject) then
            return PuiseuxSeriesObject:-DeepCopy(_passed);
        else
            return Object(PuiseuxSeriesObject, _passed);
        end if;
    end proc;

    # To get the internal power series of a Puiseux series.
    export GetPowerSeries := proc(obj :: PuiseuxSeriesObject, $)
            return obj:-GetPowerSeries(obj);
    end proc;

    # Get order of the Puiseux series variables.
    export GetPuiseuxSeriesOrder := proc(obj :: PuiseuxSeriesObject, $)
        return obj:-GetPuiseuxSeriesOrder(obj); 
    end proc;

    # Get order of the power series variables.
    export GetPowerSeriesOrder := proc(obj :: PuiseuxSeriesObject, $)
        return obj:-GetPowerSeriesOrder(obj); 
    end proc;

    # Get order of the Puiseux series.
    export GetOrder := proc(obj :: PuiseuxSeriesObject, bnd :: {nonnegint, identical(infinity)}, $)
        return obj:-GetOrder(_passed); 
    end proc;

    # Get the monomial that multiplies the PuSO.
    export GetMonomial := proc(obj :: PuiseuxSeriesObject, $)
        return obj:-GetMonomial(obj); 
    end proc;

    # Get the rays.
    export GetRays := proc(obj :: PuiseuxSeriesObject, $)
        return obj:-GetRays(obj); 
    end proc;

    # To set a nonzero_pso_bound or a nonzero_pso_bound_static.
    export SetNonzeroPowerSeriesDegreeBound := proc(bound::nonnegint,
                                            obj :: PuiseuxSeriesObject := NULL,
                                            $)
        if obj=NULL then
            return obj:-SetNonzeroPowerSeriesDegreeBound(PuiseuxSeriesObject,
                                                         bound, 'instance' = false);
        else
            return obj:-SetNonzeroPowerSeriesDegreeBound(obj, bound); 
        end if;
    end proc;

    # To set a smallest_pso_bound or a smallest_pso_bound_static.
    export SetSmallestTermDegreeBound := proc(bound::nonnegint,
                                              obj :: PuiseuxSeriesObject := NULL,
                                              $)
        if obj=NULL then
            return obj:-SetSmallestTermDegreeBound(PuiseuxSeriesObject,
                                                   bound, 'instance' = false);
        else
            return obj:-SetSmallestTermDegreeBound(obj, bound); 
        end if;
    end proc;

    # To set a Puiseux_bound or a Puiseux_bound_static.
    export SetPuiseuxBound := proc(bound::nonnegint,
                                   obj :: UnivariatePolynomialOverPowerSeriesObject := NULL,
                                   $)
        if obj=NULL then
            return obj:-SetPuiseuxBound(UnivariatePolynomialOverPowerSeriesObject,
                                        bound, 'instance' = false);
        else
            return obj:-SetPuiseuxBound(obj, bound); 
        end if;
    end proc;

    # Get the change of variables applied to the internal ps.
    export ChangeOfVariables := proc(obj :: PuiseuxSeriesObject, $)
        return obj:-ChangeOfVariables(obj); 
    end proc;

$include "MultivariatePowerSeries/PowerSeries/src/PowerSeries.mm"
$include "MultivariatePowerSeries/UPoPS/src/UPoPS.mm"
$include "MultivariatePowerSeries/PuiseuxSeries/src/PuiseuxSeries.mm"
$include "MultivariatePowerSeries/LaurentSeries/src/LaurentSeries.mm"
 
export
    LaurentSeries := proc()
        return Object(LaurentSeriesObject, _passed);
        
# OdeManipulator sub-package. For now, all the commands in this 
# module are hidden, i.e., they are not exported.
$include "MultivariatePowerSeries/OdeManipulator/src/OdeManipulator.mpl"

$undef PSDS_ENTRYTYPE 
$undef PSDS_TYPE 
$undef UPOPSDS_ENTRYTYPE
$undef UPOPSDS_TYPE

end module:

protect( 'MultivariatePowerSeries' ):
#savelib( 'MultivariatePowerSeries' ):
 
