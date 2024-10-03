######################
# ## Basic Routines  
######################

# DeepCopy
# to deep copy of the input object
export 
    DeepCopy ::static := proc(self :: PowerSeriesObject, $)
        return Object(PowerSeriesObject, copy(self:-hpoly), self:-deg, self:-gen, self:-vars,
                      self:-ancestors, self:-algexpr, self:-dstyle);
    end proc;

# SetHomogeneousPart
# to initialize idx-th Homogeneous Part of self with value:: {polynom, algebraic}.
# Note this method shouldn't be documented and used internally in UPoPS/src/weierstrass_preparation.mm
export 
    SetHomogeneousPart ::static := proc(self :: PowerSeriesObject, 
                                        idx :: nonnegint,
                                        value :: {polynom, algebraic},
                                        $)
        local dummy;
        if idx > self:-deg then
            self:-hpoly(idx+1) := 0;
            self:-deg := idx;
            self:-hpoly[idx] := AUTO_EXPAND(dummy, value);
        else
            self:-hpoly[idx] := AUTO_EXPAND(dummy, value);
        end if;
    end proc;

export
    UpdatePrecision :: static := proc(self :: PowerSeriesObject,
                                      deg :: nonnegint,
                                      $)
        HOMOGENEOUS_PART(self, deg);
        return self;
    end proc;

# GetCoefficient
# to return the coefficient of monomial mon in hpoly.
# mon : a monomial 
export 
    GetCoefficient ::static := proc(self :: PowerSeriesObject,
                                      mon :: polynom,
                                      $)
        if mon = 1 then 
            return HOMOGENEOUS_PART(self, 0);
        end if;
        local facts := convert(mon, 'list', ':-`*`');
        local pairs := map(convert, facts, 'list', ':-`^`');
        if not type(pairs, 'list'(['name', 'integer'])) then
            error "invalid input: expected monomial, but found %1", monomial;
        end if;
        local j;
        local d := add(j[2], j in pairs);
        self:-ensure_degree(self, d);
        #find the term:
        local c := self:-hpoly[d];
        for local i from 1 to numelems(pairs) do
            c := coeff(c, pairs[i][1], pairs[i][2]);
        end do;
        return c;
    end proc;

# GetPrecision
# to return the upperbound of poly:-hpoly 
export 
    GetPrecision ::static := proc(self :: PowerSeriesObject, $) 
        return self:-deg;
    end proc;

# GetAncestor
# to get the ancestor named aname from self (self:-ancestors[aname])
# Note this method shouldn't be documented and used internally only in UPoPS/src/weierstrass_preparation.mm
export 
    GetAncestor ::static := proc(self :: PowerSeriesObject,      
                                aname,
                                $)
        if aname in {exports(self:-ancestors)} then
            return self:-ancestors[aname];
        else
            error "ancestor name %1 doesn't exist", aname;
        end if;
    end proc;

# GetAnalyticExpression
# to get the analytic expression of the input object (self:-algexpr)
export 
    GetAnalyticExpression ::static := proc(self :: PowerSeriesObject, $) 
        return self:-algexpr;
    end proc; 

# HomogenousPart
# to compute the PowerSeriesObject to degree d (if d < self:-deg) and return self:-hpoly[d]  
# d : the index (totaldegree) of required homogeous part of the input object
export 
    HomogeneousPart ::static := proc(self :: PowerSeriesObject,
                                    d :: nonnegint,
                                    $)
        local dummy;
        return AUTO_EXPAND(dummy, HOMOGENEOUS_PART(self, d));
    end proc;

# shouldn't be documented! 
export 
    HomogeneousPart_local ::static := proc(self :: PowerSeriesObject,
                                    d :: nonnegint,
                                    $)
        return HOMOGENEOUS_PART(self, d);
    end proc;

# Zero
# to create Zero PowerSeriesObject
export 
    Zero ::static := proc()
        local zero_gen := proc(self :: PowerSeriesObject, 
                               d :: nonnegint,
                               $)
            return 0;
        end proc; 
        return Object(PowerSeriesObject, Array(0 .. 0, [0]), 0, zero_gen, {}, 0);
    end proc;

# One
# to create One PowerSeriesObject
export 
    One ::static := proc()
        local one_gen := proc(self :: PowerSeriesObject, 
                            d :: nonnegint,
                            $)
            return ifelse(d = 0, 1, 0);
        end proc;
        return Object(PowerSeriesObject, Array(0 .. 0, [1]), 0, one_gen, {}, 1);
    end proc;

# Identity
# to create Identity PowerSeriesObject
export
    Identity :: static := proc(x :: name, $)
        local iden_gen := proc(self :: PowerSeriesObject, 
                            d :: nonnegint,
                            $)
            return ifelse(d = 1, self:-hpoly[1], 0);
        end proc;

        return Object(PowerSeriesObject, Array(0 .. 1, [0, x]), 1, iden_gen, {x}, x);
    end proc;

# Constant 
# to create Constant PowerSeriesObject
export 
    Constant ::static := proc(p :: COEFFICIENT_TYPE, $)
        if p = 0 then
            return Zero();
        elif p = 1 then
            return One();
        end if;
        
        local const_gen := proc(self :: PowerSeriesObject, 
                            d :: nonnegint,
                            $)
            return ifelse(d = 0, self:-hpoly[0], 0);
        end proc;

        return Object(PowerSeriesObject, Array(0 .. 0, [p]), 0, const_gen, {}, p);
    end proc;

# IsUnit
# to check if self is unit 
export 
    IsUnit ::static := proc(self :: PowerSeriesObject,
                                {evaluation::identical(floats,evala,none) := 'usenone'},
                                { tolerance :: And( numeric, nonnegative ) := Float( 1, 2 - Digits ) },
                                { guarddigits :: nonnegint := 3 }, $)
        local dummy, alpha;
        local my_hom_part := HOMOGENEOUS_PART(self, 0);
        local expand_dummy := AUTO_EXPAND(dummy,my_hom_part);

        if evaluation=':-floats' then  
            alpha := evalf[Digits+guarddigits](expand_dummy);
            alpha := fnormal( alpha, Digits, tolerance );
        elif evaluation=':-evala' then 
            alpha := evala(simplify(expand_dummy));
        else
            alpha := expand_dummy
        end if;
        
        return evalb(alpha <> 0);

    end proc;

# ApproximatelyZero
# to check approximately (up to a specific precision) if self is zero
# if deg = -2, then this algorithm will compare ONLY generated hpolys up to self:-deg with zero polynomial,
# elif deg = -1, returns always true, else, compare self:-hpoly[i] and 0 for i = 0, ..., deg.
export 
    ApproximatelyZero ::static := proc(self :: PowerSeriesObject,
                            deg :: integer := -2,
                            {force :: truefalse := false},
                            $)
        # if the algexpr equals 0
        if (not force) and self:-algexpr <> undefined and self:-algexpr = 0 then 
            return true;
        end if;
        #
        local d := ifelse(deg < -1, self:-deg, deg);
        for local i from 0 to d do
            self:-ensure_degree(self, i);
            if self:-hpoly[i] <> 0 then
                return false;
            end if;
        end do;
        return true;
    end proc;

# ApproximatelyEqual
# to check approximately the equality of two PowerSeriesObject
# deg is the maximum precision; base on that self and other will be computed and compared 
# if deg = -2, then this algorithm will compare the first max(self:-deg, other:-deg) homogeneous parts,
# elif deg = -1, returns always true, else, compare self and other up to deg.
export 
    ApproximatelyEqual ::static := proc(self :: PowerSeriesObject, 
                            other :: PowerSeriesObject,
                            deg :: integer := -2,
                            {force :: truefalse := false},
                            $)
        # if the algexprs equal 
        if (not force) and self:-algexpr <> undefined and self:-algexpr = other:-algexpr then 
            return true;
        end if;
        # 
        local d := ifelse(deg < -1, max(self:-deg, other:-deg), deg);
        for local i from 0 to d do
            self:-ensure_degree(self, i);
            other:-ensure_degree(other, i);
            if self:-hpoly[i] <> other:-hpoly[i] then
                return false;
            end if;
        end do;
        return true; 
    end proc;

# FromProcedure
# to create a PowerSeriesObject from a procedure or an appliable_module
# p : a procedure or appliable_nodule with one integer/nonnegint input (e.g. d -> x^d)
# expr : if defined, then set the algexpr 
# validity_check : if true, then check the validity of the results of p 
# do_expand : if true, then expand the results of p 
export 
    FromProcedure ::static := proc(p :: {procedure, appliable_module}, 
                                   expr :: algebraic := undefined,
                                   vars :: {set(name), identical(undefined)} := undefined,
                                   {check :: truefalse := true},
                                   {expand :: truefalse := true},
                                   $)
        local ps, the_vars;
        if expr <> undefined then
            if type(vars, ':-set') then
                local missing := indets(expr, ':-name') minus vars;
                missing := select(x -> depends(expr, x), missing);
                if missing <> {} then
                    error "invalid input: the given analytic expression, %1, contains variables %2 which were "
                    "not specified", expr, missing;
                end if;
                the_vars := vars;
            else
                the_vars := indets(expr, ':-name');
            end if;
            ps := Object(PowerSeriesObject, Array(0 .. 0), 0, proc_gen, the_vars, ["P" = p, "V" = check, "E" = expand], expr);
        elif vars = undefined then
            error "invalid input: if you don't specify the analytic expression for procedure input, "
            "you must specify the set of variables";
        else
            ps := Object(PowerSeriesObject, Array(0 .. 0), 0, proc_gen, vars, ["P" = p, "V" = check, "E" = expand]);
        end if;
        ps:-hpoly[0] := proc_gen(ps, 0); 
        return ps;
    end proc;

# FromPolynomial
# to create a PowerSeriesObject from a polynomial
# p : a polynomial 
export 
    FromPolynomial ::static := proc(p :: polynom, $)
        local new_poly;
        if not type(p, ':-polynom'(':-complex'(':-numeric'))) then
            new_poly := Algebraic:-Expand(p);
        else 
            new_poly := expand(p);
        end if;
        if type(new_poly, 'COEFFICIENT_TYPE') then
            return Constant(new_poly);
        end if;
        local homo_poly_array := convert_hpoly_from_poly(new_poly);

        return Object(PowerSeriesObject, homo_poly_array, degree(new_poly), poly_gen, 
                      indets( new_poly, ':-And'(':-name',':-Not'(':-constant'))), p);
    end proc;

# Truncate
# to convert to a polynomial from a PowerSeriesObject up to degree (precision) deg
# deg : if deg < 0, convert up to degree self:-deg (default), else convert up to degree deg
export 
    Truncate ::static := proc(self :: PowerSeriesObject, 
                                deg :: integer := -1,
                                $)
        if deg < 0 then
            return add(self:-hpoly); 
        end if;
        self:-ensure_degree(self, deg);
        return add(self:-hpoly[0 .. deg]);
    end proc;

# GeometricSeries
# to create a geometric series
# vars : a list of variables (names)  
export 
    GeometricSeries ::static := proc(vars :: {name, nonemptylist(name)}, $)
        if type(vars, ':-name') then
            return thisproc([vars]);
        end if;
        if numelems(vars) <> numelems(convert(vars, ':-set')) then
            error "invalid input: the variables have to be distinct";
        end if;
        local geo_gen := proc(self :: PowerSeriesObject, 
                            d :: nonnegint,
                            $) 
            local i;
            local s := add(vars[i], i = 1 .. nops(vars));
            return expand(s^d);
        end proc;
        local v;
        return Object(PowerSeriesObject, Array(0 .. 1, [1, add(vars)]), 1, geo_gen, {op(vars)}, 1/(1 - add(v, v in vars))); 
    end proc;

# SumOfAllMonomials
# to create a sum_of_all_monomials
# vars : a list of variables (names)
export
    SumOfAllMonomials :: static := proc(vars :: {name, nonemptylist(name)}, $)
        if type(vars, ':-name') then
            return thisproc([vars]);
        end if;
        if numelems(vars) <> numelems(convert(vars, ':-set')) then
            error "invalid input: the variables have to be distinct";
        end if;
        local resulting_object;
        local sam_gen := proc(self :: PowerSeriesObject, 
                              d :: nonnegint, 
                              $)
            if d = 0 then
                return 1;
            end if;
            local previous_degree := HomogeneousPart(resulting_object, d-1);
            local terms := Array(1 .. numelems(vars));
            for local i from numelems(vars) to 1 by -1 do
                terms[i] := vars[i] * previous_degree;
                previous_degree := subs(vars[i] = 0, previous_degree);
            end do;
            return expand(add(terms));
        end proc;
        local v;
        resulting_object := Object(PowerSeriesObject, Array(0..1, [1, add(vars)]), 1, sam_gen, {op(vars)}, 1/mul(1-v, v in vars)); 
        return resulting_object;
    end proc;

export
    Variables ::static := proc(self :: PowerSeriesObject, $)
        return self:-vars;
    end proc;

local
monomial_coefficient_gen ::static := proc(self :: PowerSeriesObject,
                                          d :: nonnegint,
                                          $)
local hom_part := HOMOGENEOUS_PART(self:-ancestors:-A, d + self:-ancestors:-D);
    return coeff(hom_part, self:-ancestors:-X, self:-ancestors:-D);
end proc;

# The coefficient of x^d in self, as a power series. For use in UPoPS.
export
MonomialCoefficient ::static := proc(self :: PowerSeriesObject,
                                     x :: name,
                                     d :: nonnegint,
                                     $)
    if not x in self:-vars then
        if d = 0 then
            return self;
        else
            return Zero();
        end if;
    elif self:-vars = {x} then
        return Constant(eval(HOMOGENEOUS_PART(self, d), x = 1));
    end if;

local algexpr;
    if self:-algexpr = undefined then
        algexpr := undefined;
    elif type(self:-algexpr, ':-polynom'(':-anything', x)) then
        if d > degree(self:-algexpr, x) then
            return Zero();
        end if;
        algexpr := coeff(self:-algexpr, x, d);
    else
        algexpr := eval(diff(self:-algexpr, x $ d), x = 0) / d!;
    end if;

local first_entry := coeff(HOMOGENEOUS_PART(self, d), x, d);
    return Object(PowerSeriesObject, Array(0..0, [first_entry]), 0, monomial_coefficient_gen,
                  self:-vars minus {x}, ["A" = self, "D" = d, "X" = x], algexpr);
end proc;

export
ConvertToUPoPS ::static := proc(self :: PowerSeriesObject,
                                x :: name,
                                $)
    if not member(x, self:-vars) then
        return Object(UnivariatePolynomialOverPowerSeriesObject, [self], x);
    end if;
    
    if self:-algexpr = undefined or not type(self:-algexpr, ':-polynom'(':-anything', x)) then
        error "attempted to convert a power series involving %1 to a univariate polynomial over "
        "power series in %1, but it is not known to be polynomial in %1", x;
    end if;
    
local i, d := degree(self:-algexpr, x);
    return Object(UnivariatePolynomialOverPowerSeriesObject,
                  [seq(self:-MonomialCoefficient(self, x, i), i = 0 .. d)], x);
end proc;
