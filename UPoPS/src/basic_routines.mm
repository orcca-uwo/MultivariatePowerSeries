#######################################################
## # UPoPS Basic Routines
#######################################################

# DeepCopy
# to deep copy of the input object
export 
    DeepCopy ::static := proc(self :: UnivariatePolynomialOverPowerSeriesObject, $)
        local new_upoly := map(function:-DeepCopy, self:-upoly);
        local upops := Object(UnivariatePolynomialOverPowerSeriesObject, new_upoly, self:-vname);

        upops:-Puiseux_bound := self:-Puiseux_bound;
        upops:-Hensel_bound := self:-Hensel_bound;

        return upops;
    end proc;

# GetCoefficient
# to get the d-th coefficient of the input object
# d : nonnegint
export 
    GetCoefficient ::static := proc(self :: UnivariatePolynomialOverPowerSeriesObject, 
                                    d :: nonnegint,
                                    $)
        if DEGREE(self) < d then 
            error "invalid input: degree out of range";
        else 
            return self:-upoly[d];
        end if;
    end proc;

# GetAnalyticExpression
# to get the analytic expression of the input object
export 
    GetAnalyticExpression ::static := proc(self :: UnivariatePolynomialOverPowerSeriesObject, $)
        local i, lt;
        local algexprs := [seq(function:-GetAnalyticExpression(lt), lt in self:-upoly)];
        if member(undefined, algexprs) or self:-vname = undefined then 
            return 'undefined';
        else 
            return add(algexprs[i+1] * (self:-vname)^i, i = 0 .. DEGREE(self));
        end if;
    end proc;

# Degree
# to get the degree of the input object (DEGREE(self))
export 
    Degree ::static := proc(self :: UnivariatePolynomialOverPowerSeriesObject, $)
        return DEGREE(self);
    end proc;

# MainVariable
# to get the main variable of self, self:-vname if it's definded 
export 
    MainVariable ::static := proc(self :: UnivariatePolynomialOverPowerSeriesObject, $)
        if self:-vname <> undefined then 
            return self:-vname;
        else 
            return 'undefined';
        end if;
    end proc;

# UpdatePrecision
# to ensure the coefficients of self :: UPoPS are computed up tp degree 'prec'. 
# prec : nonnegint 
export 
    UpdatePrecision ::static := proc(self :: UnivariatePolynomialOverPowerSeriesObject, 
                                    prec :: nonnegint, 
                                    $)
    for local i from 0 to DEGREE(self) do
        if self:-IsPuSOUPoP(self)=true then
            PowerSeriesObject:-HomogeneousPart_local(self:-upoly[i]:-GetPowerSeries(self:-upoly[i]), 
                                                        prec);
        else 
            PowerSeriesObject:-HomogeneousPart_local(self:-upoly[i], prec);
        end if;
    end do;
    return self;
end proc;

# FromPolynomial
# to get a UPoPS from a polynomial w.r.t the main variable x::name
# x : main variable 
export 
    FromPolynomial ::static := proc(p :: polynom, 
                            x :: name,
                            $)
        local new_p, dummy;
        if type(p, 'COEFFICIENT_TYPE') then
            return UnivariatePolynomialOverPowerSeriesObject:-Constant(p);
        else
            new_p := AUTO_EXPAND(dummy, p); 
        end if;
        # local vars := indets(new_p);
        # if not member(x, vars) then
        #     error "invalid input: variable %1 must be in the variable list, %2", x, vars;
        # end if;
        local d := degree(new_p, x);
        local cfs := PolynomialTools:-CoefficientList(new_p, x);
        local upoly := Array(0 .. d, i -> PowerSeriesObject:-FromPolynomial(cfs[i+1]));
        
        return Object(UnivariatePolynomialOverPowerSeriesObject, upoly, x);
    end proc;

# ConvertToPowerSeries 
# to convert a UPoPS to a PowerSeriesObject
export 
    ConvertToPowerSeries ::static := proc(self :: UnivariatePolynomialOverPowerSeriesObject, $)
        if self:-IsPuSOUPoP(self)=true then
            error "invalid input for the ConvertToPowerSeries function %1 must "
                    "have power series coefficients.", up;
        end if;

        local i;
        local coefs := [seq(self:-vname^i, i = 0 .. DEGREE(self))];
        
        return PowerSeriesObject:-NaryAdd(seq(self:-upoly), ':-coefficients' = coefs);
    end proc;

# truncate a UPoPS to polynomial e.g. UPoPSObject:-ToPolynomial 
# depends on the mode,
# mode = totaldegree : truncate self up to totaldegree = deg 
# mode = powerseriesdegree : update the coefficients of self up to precision deg and return the converted polynomial
# Note mode = powerseriesdegree is default
export 
    Truncate ::static := proc(self :: UnivariatePolynomialOverPowerSeriesObject, 
                            deg :: integer := -1,
                            {mode :: identical(totaldegree, powerseriesdegree) := powerseriesdegree},
                            $)
        local dummy, poly := 0;
        local x := ifelse(self:-vname <> undefined, self:-vname, '_Z');
        if mode = ':-totaldegree' then
            if deg < -1 then
                error "invalid input: cannot truncate at degree %1 in totaldegree mode", deg;
            end if;
            for local i from 0 to min(DEGREE(self), deg) do 
                if type(self:-upoly[i], PowerSeriesObject) then
                    poly += AUTO_EXPAND(dummy, function:-Truncate(self:-upoly[i], deg - i ) * x^i );
                else 
                    poly += expand(function:-Truncate(self:-upoly[i], deg) * x^i); 
                end if;
            end do;
        else 
            for local i from 0 to DEGREE(self) do
                if type(self:-upoly[i], PowerSeriesObject) then
                    poly += AUTO_EXPAND(dummy, function:-Truncate(self:-upoly[i], deg) * x^i); 
                else 
                    poly += expand(function:-Truncate(self:-upoly[i], deg) * x^i); 
                end if;
            end do;
        end if;
        return collect(poly, x);
    end proc;

# Zero
# to create a zero UPoPS
# Note: UPoPS:-Zero is 0 and not an empty UPoPS w/ Array(0 .. -1), we use this feature in basic_arithmetic and WP. 
export 
    Zero ::static := proc()
        return Object(UnivariatePolynomialOverPowerSeriesObject, Array(0 .. 0, [PowerSeriesObject:-Zero()])); 
    end proc;

# One 
# to create a one UPoPS
export 
    One ::static := proc()
        return Object(UnivariatePolynomialOverPowerSeriesObject, Array(0 .. 0, [PowerSeriesObject:-One()]));
    end proc;

# Constant 
# to create a constant UPoPS
export 
    Constant ::static := proc(p :: COEFFICIENT_TYPE, $)
        if p = 0 then
            return Zero();
        elif p = 1 then
            return One();
        else
            return Object(UnivariatePolynomialOverPowerSeriesObject, Array(0 .. 0, [PowerSeriesObject:-Constant(p)]));
        end if;
    end proc;


# ApproximatelyZero 
# see PowerSeriesObject:-ApproximatelyZero for details 
export 
    ApproximatelyZero ::static := proc(self :: UnivariatePolynomialOverPowerSeriesObject,
                            deg :: integer := -2,
                            {force :: truefalse := false},
                            $)    
    local sw := self:-IsPuSOUPoP(self);
                            
    for local i from 0  to upperbound(self:-upoly) do
        if sw = false then
            if not PowerSeriesObject:-ApproximatelyZero(self:-upoly[i], deg, _options['force']) then
                return false;
            end if;
        else 
            if not PuiseuxSeriesObject:-ApproximatelyZero(self:-upoly[i], deg) then
                return false;
            end if;
        end if;
    end do;
    return true;
end proc;

# IsUnit
# to check if the self:-upoly[0] is a unit
export 
    IsUnit ::static := proc(self :: UnivariatePolynomialOverPowerSeriesObject, $)
        if self:-IsPuSOUPoP(self)= true then 
            error "invalid input for the IsUnit function %1 must "
                  "have power series coefficients.", self;
        else
            return PowerSeriesObject:-IsUnit(self:-upoly[0]);
        end if;
    end proc;

# ApproximatelyEqual
# see PowerSeriesObject:-ApproximatelyEqual for details
export 
    ApproximatelyEqual ::static := proc(_self :: UnivariatePolynomialOverPowerSeriesObject, 
                            _other :: UnivariatePolynomialOverPowerSeriesObject,
                            deg :: integer := -2,
                            {force :: truefalse := false},
                            $)
    local self, other;
    local sw := false;
    if _self:-IsPuSOUPoP(_self)= true or _other:-IsPuSOUPoP(_other)= true then        
        self := _self:-ConvertToPuSOUPoP(_self);
        other := _other:-ConvertToPuSOUPoP(_other);
        sw := true;
    else 
        self := _self;
        other := _other;
    end if;

    if deg = -1 then 
        return true;
    elif DEGREE(self) <> DEGREE(other) or upperbound(self:-upoly) <> upperbound(other:-upoly) then
        return false;
    end if;
    if upperbound(self:-upoly) < 0 then return true; end if;
    for local i from 0 to DEGREE(self) do 
        if sw = false then
            if not PowerSeriesObject:-ApproximatelyEqual(self:-upoly[i], other:-upoly[i], deg, _options['force']) then
                return false;
            end if;
        else 
            if not PuiseuxSeriesObject:-ApproximatelyEqual(self:-upoly[i], other:-upoly[i], deg) then
                return false;
            end if;
        end if;
    end do;
    return true;
end proc;

# a local function to return the common vname of self and other 
# mainly used in basic_arithmetic methods 
# shouldn't be documented 
local 
    get_common_vname :: static := proc(objs :: seq(UnivariatePolynomialOverPowerSeriesObject),
                                       $)
    local o, vnames := {seq(o:-vname, o in [objs])} minus {undefined};
        if numelems(vnames) = 0 then
            return undefined;
        elif numelems(vnames) = 1 then
            return vnames[1];
        else
            error "incompatible inputs: expect UnivariatePolynomialOverPowerSeries with the "
            "same main variable, but received %1 <> %2", vnames[1], vnames[2];
        end if;
    end proc;

export
    Variables ::static := proc(self :: UnivariatePolynomialOverPowerSeriesObject, $)
        if self:-vname = undefined then
            if DEGREE(self) = 0 then
                return function:-Variables(self:-upoly[0]);
            else
                error "invalid input: this univariate polynomial over power series' main variable is "
                "unset, so its set of variables is not well-defined";
            end if;
        else
            local pso;
            return `union`(seq(pso:-Variables(pso), pso in self:-upoly), {self:-vname});
        end if;
    end proc;

# To know if it is a UPOP with PuSO coefficients.
# If the coefficients of self are PuSOs then the function
# returns true. Otherwise, false is returned.
export IsPuSOUPoP ::static := proc(self :: UnivariatePolynomialOverPowerSeriesObject, $)
    return ormap(type, self:-upoly, PuiseuxSeriesObject);
end proc;

# To convert a UPoP with PSO coefficient to a UPoP with PuSO as 
# coefficients.
export ConvertToPuSOUPoP::static := proc(self :: UnivariatePolynomialOverPowerSeriesObject, 
                                            mp::list(`=`( name, {:-`*`( { name, name^rational } ), name, name^rational})) := [], $)
    if self:-IsPuSOUPoP(self) = true then
        return self;
    else
        local A := map(PuiseuxSeries, self:-upoly, mp);

        return MultivariatePowerSeries:-UnivariatePolynomialOverPowerSeries(A, self:-vname);
    end if;
end proc;

# To know if a UPoPS is monic.
export IsMonic::static := proc(self :: UnivariatePolynomialOverPowerSeriesObject, $)
    local n := DEGREE(self);
    local a_n := self:-GetCoefficient(self, n);

    local ana := a_n:-GetAnalyticExpression(a_n);

    if normal(ana)=1 then
        return true;
    else 
        return false;
    end if;
end proc;

# To make a UPoPS monic as a UPoPuS.
export MakeMonic::static := proc(_self :: UnivariatePolynomialOverPowerSeriesObject, 
                                 {returnleadingcoefficient :: {truefalse, identical(automatic)}
                                  := ':-automatic'}, $)
    local self := _self:-ConvertToPuSOUPoP(_self);
    local n := DEGREE(self);
    local a := self:-GetCoefficient(self, n);

    if self:-IsMonic(self) then 
        if returnleadingcoefficient = true then
            return [self, Object(a, 1)];
        else 
            return self;
        end if;
    end if;
    
    local a_inverse := a:-Inverse(a);
    local upoly := self:-upoly;
    upoly := map(Multiply, upoly, a_inverse);
    upoly[n] := Object(a, 1);

    local output := Object(UnivariatePolynomialOverPowerSeriesObject, upoly, self:-vname);
    
    if returnleadingcoefficient = true then
        return [output, a];
    else 
        return output;
    end if;
end proc;