######################
# ## Basic Arithmetic  
######################

# BinaryAdd
# to add two objects
export
    BinaryAdd ::static := proc(self :: PowerSeriesObject,
                                other :: PowerSeriesObject,
                                $)
        local new_deg := min(self:-deg, other:-deg);
        local new_hpoly := zip(':-`+`', UP_TO_DEGREE_ARRAY(self, new_deg), UP_TO_DEGREE_ARRAY(other, new_deg));
        return Object(PowerSeriesObject, new_hpoly, new_deg, add_gen, self:-vars union other:-vars,
                      ["A" = self, "B" = other],
                      ifelse(membertype(undefined, [self:-algexpr, other:-algexpr]), undefined,
                             self:-algexpr + other:-algexpr));
    end proc;

# NaryAdd
# given a sequence of Powerseries, polynom, or algebraic and a list of coefficients of polynom or algebraic
# return \sum_{i = 0}^{numelems([terms])} coefficients[i] * terms[i].
export 
    NaryAdd ::static := proc(terms :: seq({PowerSeriesObject, algebraic}), 
                             {coefficients :: list(algebraic) := []},
                             $)
        if numelems(coefficients) <> 0 then 
            return NaryAddWithCoefficients(_passed);
        end if;

        local lterms := [terms];
        local new_deg, new_hpoly, i, t, lt;
        if not type(lterms, 'list'(PowerSeriesObject)) then
            local nps;
            lterms, nps := selectremove(type, lterms, PowerSeriesObject);
            local nps_as_ps := FromPolynomial(add(t, t in nps));
            lterms := [op(lterms), nps_as_ps];
        end if;
        if numelems(lterms) < 2 then
            return ifelse(numelems(lterms) = 1, DeepCopy(lterms[1]), Zero());
        elif numelems(lterms) = 2 then 
            return BinaryAdd(lterms[1], lterms[2]);
        end if;
        new_deg := min(seq(t[':-deg'], t in lterms));
        new_hpoly := Array(0 .. new_deg, [seq(add(t[':-hpoly'][i], t in lterms) , i = 0 .. new_deg )]);
        local algexprs := [seq(lt:-algexpr, lt in lterms)];
        local expr := ifelse(membertype(undefined, algexprs), undefined, add(algexprs));
        local vars := `union`(seq(lt:-vars, lt in lterms));
        return Object(PowerSeriesObject, new_hpoly, new_deg, nary_add_gen, vars, ["T" = lterms], expr);
    end proc;

# local function called in PowerSeriesObject:-`+` 
# given a sequence of Powerseries, polynom, or algebraic and a list of coefficients of polynom or algebraic
# return \sum_{i = 0}^{numelems([terms])} coefficients[i] * terms[i].
# Note Shoudn't be documented 
local 
    NaryAddWithCoefficients ::static := proc(terms :: seq({PowerSeriesObject, algebraic}), 
                                             {coefficients :: list(algebraic) := []},
                                             $)
        local lterms := [terms];
        local n := numelems(lterms);
        if n <> numelems(coefficients) then
            error "invalid input: coefficients, when specified, should have as many members "
            "as there are arguments, but received %1 arguments with coefficients = %2", n, coefficients;
        end if;

        local pairs := zip(`[]`, lterms, coefficients);
        pairs := remove(type, pairs, [':-anything', 0]);

        local nps, dummy;
        pairs, nps := selectremove(type, pairs, [PowerSeriesObject, anything]);
        if nps <> [] then
            local t, nps_as_ps := FromPolynomial(add(AUTO_EXPAND(dummy, t[1] * t[2]), t in nps));
            pairs := [op(pairs), [nps_as_ps, 1]];
        end if;

        pairs := map(pair -> (if type(pair[2], {':-numeric', ':-algnum', ':-algnum'^':-fraction'}) then
                                  pair;
                              else
                                  [BinaryMultiply(pair[1], FromPolynomial(pair[2])), 1];
                              end if), pairs);

        local p, new_deg := min(seq(p[1][':-deg'], p in pairs));
        local i, new_hpoly := Array(0 .. new_deg, [
            seq(add(AUTO_EXPAND(dummy, p[1][':-hpoly'][i] * p[2]), p in pairs), i = 0 .. new_deg)]);

        local lt, algexprs := [seq(lt[1][':-algexpr'] * lt[2], lt in pairs)];
        local expr;
        if membertype(undefined, algexprs) then
            expr := undefined;
        else
            expr := add(algexprs);
        end if;

        local vars := `union`(seq(lt:-vars, lt in lterms));

        return Object(PowerSeriesObject, new_hpoly, new_deg, nary_add_with_coeffs_gen, vars, ["P" = pairs], expr);
    end proc;

# Powerseries:-`+` n-ary operator 
export 
    `+` ::static := proc(terms :: seq({PowerSeriesObject, algebraic}), $)
            return NaryAdd(terms);
    end proc;

# BinarySub
# to subtract two objects
export
    BinarySub ::static := proc(self :: PowerSeriesObject,
                                other :: PowerSeriesObject,
                                $)
        local new_deg := min(self:-deg, other:-deg);
        local new_hpoly := zip(':-`-`', UP_TO_DEGREE_ARRAY(self, new_deg), UP_TO_DEGREE_ARRAY(other, new_deg));
        return Object(PowerSeriesObject, new_hpoly, new_deg, sub_gen, self:-vars union other:-vars,
                      ["A" = self, "B" = other], ifelse(membertype(undefined, [self:-algexpr, other:-algexpr]),
                                                        undefined, self:-algexpr - other:-algexpr));
    end proc;

# Negate
# to negate the input object
export 
    Negate ::static := proc(self :: PowerSeriesObject, $)
        return Object(PowerSeriesObject, - self:-hpoly, self:-deg, neg_gen, self:-vars, ["A" = self],
                      ifelse(type(self:-algexpr, undefined), undefined, - self:-algexpr));
    end proc;

# BinaryMultiply
# to multiply two objects 
export 
    BinaryMultiply ::static := proc(self :: PowerSeriesObject,
                                    other :: PowerSeriesObject,
                                    $)
        local dummy;
        local new_deg := min(self:-deg, other:-deg);
        local a := self:-hpoly;
        local b := other:-hpoly;
        local j;
        local new_hpoly := Array(0 .. new_deg, i -> (AUTO_EXPAND(dummy, add(a[j] * b[i-j], j = 0 .. i))));
        return Object(PowerSeriesObject, new_hpoly, new_deg, mul_gen, self:-vars union other:-vars,
                      ["A" = self, "B" = other],
                      ifelse(membertype(undefined, [self:-algexpr, other:-algexpr]), undefined,
                             self:-algexpr * other:-algexpr));
    end proc;

# a local function used in NaryMultiply when size=constant
local
    multiply_binary_split :: static := proc(lterms :: nonemptylist(PowerSeriesObject), $)
        if numelems(lterms) = 1 then
            return lterms[1];
        else
            local k := trunc(numelems(lterms)/2);
            return BinaryMultiply(thisproc(lterms[1..k]), thisproc(lterms[k+1..]));
        end if;
    end proc;

# NaryMultiply
# to multiply a sequence of power series, polynomials, and algebraics
# size : if the degree pattern is known, then use one of these three different modes 
#        {constant, increasing, decreasing} or leave it as default 
export
    NaryMultiply ::static := proc(factors :: seq({PowerSeriesObject, algebraic}), 
                                {size :: identical(constant, increasing, decreasing) := constant},
                                $)
        local lfactors := [factors];
        local n := numelems(lfactors);
        local t;
        if not type(lfactors, 'list'(PowerSeriesObject)) then
            local nps;
            lfactors, nps := selectremove(type, lfactors, PowerSeriesObject);
            local nps_as_ps := FromPolynomial(mul(t, t in nps));
            lfactors := [op(lfactors), nps_as_ps];
        end if;
        n := numelems(lfactors);
        if n < 2 then
            return ifelse(n = 1, DeepCopy(lfactors[1]), One());
        elif size = ':-constant' then 
            # recursively spliting 
            return multiply_binary_split(lfactors);
        elif size = ':-increasing' then 
            # f(f(f(a, b), c), d)
            return foldl(BinaryMultiply, seq(lfactors));
        else  
            # f(a, f(b, f(c, d)))
            return foldr(BinaryMultiply, lfactors[-1], op(.. -2,lfactors));
        end if;
    end proc;

# Powerseries:-`*` n-ary operator 
export 
    `*` ::static := proc(factors :: seq({PowerSeriesObject, algebraic}), $)
        return NaryMultiply(factors);
    end proc;

# BinaryExactQuotient
# to compute the exact quotient of two objects 
export 
    BinaryExactQuotient ::static := proc(self :: PowerSeriesObject, 
                                    other :: PowerSeriesObject,
                                    $)
        if HOMOGENEOUS_PART(other, 0) = 0 then 
            error "invalid input: not invertible!"; 
        end if;
        local dummy;
        local hpoly := Array(0 .. 0, [AUTO_EXPAND(dummy, normal(HOMOGENEOUS_PART(self, 0) / HOMOGENEOUS_PART(other, 0)))]);
        return Object(PowerSeriesObject, hpoly, 0, beq_gen, self:-vars union other:-vars,
                      ["A" = self, "B" = other],
                      ifelse(membertype(undefined, [self:-algexpr, other:-algexpr]), undefined,
                             self:-algexpr / other:-algexpr));
    end proc;

# Inverse
# to compute inverse, using BinaryExactQuotient 
export 
    Inverse ::static := proc(self :: PowerSeriesObject, $)
        if IsUnit(self) then 
            return BinaryExactQuotient(PowerSeriesObject:-One(), self);
        else 
            error "not invertible";
        end if;
    end proc;

# Exponentiate
# to compute self^n
# n :: negative/positive integer
export 
    Exponentiate ::static := proc(self :: PowerSeriesObject, 
                                n :: integer,
                                $)
        if n = 0 then 
            return One(); 
        elif n = 1 then 
            return DeepCopy(self);
        end if;        
        local p := ifelse(sign(n) > 0, self, Inverse(self));
        local m := abs(n);
        local res := Array(0 .. -1);
        while m > 0 do 
            if irem(m,2) = 1 then 
                res ,= p;
            end if;
            p := BinaryMultiply(p, p);
            m := iquo(m,2);
        end do;
        return NaryMultiply(seq(res));
    end proc;

# PowerSeriesObject:-`^` 
export
    `^` ::static := proc(x, y, $)
        if type(x, PowerSeriesObject) and type(y, integer) then 
            return x:-Exponentiate(x, y);
        else 
            error "invalid input: this form of exponentiation is not implemented";
        end if;
    end proc;
