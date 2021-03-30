##########################################################
# ## UPoPS Basic Arithmetic  
##########################################################

# BinaryAdd
# to add two objects
export 
    BinaryAdd ::static := proc(self :: UnivariatePolynomialOverPowerSeriesObject,
                                other :: UnivariatePolynomialOverPowerSeriesObject,
                                $)
        local x := get_common_vname(self, other);
        local u := max(DEGREE(self), DEGREE(other));
        local l := min(DEGREE(self), DEGREE(other));
        local apoly := self:-upoly;
        local bpoly := other:-upoly;
        local new_upoly := Array(0 .. l, i -> PowerSeriesObject:-BinaryAdd(apoly[i], bpoly[i]));
        if u > l then
            new_upoly(u+1) := 0; # reallocate Array 
            local tmp_upoly := ifelse(DEGREE(self) > DEGREE(other), apoly, bpoly);
            for local i from l+1 to u do
                new_upoly[i] := PowerSeriesObject:-DeepCopy(tmp_upoly[i]);
            end do;
        end if;
        return Object(UnivariatePolynomialOverPowerSeriesObject, new_upoly, x);
    end proc;

# NaryAdd 
# to add a sequence of UPoPS objects
export 
    NaryAdd ::static := proc(terms :: seq({UnivariatePolynomialOverPowerSeriesObject, PowerSeriesObject, algebraic}), $)
    local lterms := [terms];
    local vname;
        
        if not type(lterms, 'list'(UnivariatePolynomialOverPowerSeriesObject)) then
            local upops, rest;
            upops, rest := selectremove(type, lterms, UnivariatePolynomialOverPowerSeriesObject);
            vname := get_common_vname(op(upops));
            
            local pso := PowerSeriesObject:-NaryAdd(op(rest));
            lterms := [op(upops), pso:-ConvertToUPoPS(pso, vname)];
        else
            vname := get_common_vname(op(lterms));
        end if;
        
    local n := numelems(lterms);
        if n < 2 then
            return ifelse(n = 1, DeepCopy(lterms[1]), UnivariatePolynomialOverPowerSeriesObject:-Zero());
        end if;
        
    local u, d := max(seq(DEGREE(u), u in lterms));
    local upoly := Array(0 .. d);
        for local i from 0 to d do
            lterms := select(u -> u:-Degree(u) >= i, lterms);
            upoly[i] := PowerSeriesObject:-NaryAdd(seq(u:-upoly[i], u in lterms));
        end do;

        return Object(UnivariatePolynomialOverPowerSeriesObject, upoly, vname);
    end proc;

export
    `+` ::static := proc(terms :: seq({UnivariatePolynomialOverPowerSeriesObject, PowerSeriesObject, algebraic}), $)
        return NaryAdd(terms);
    end proc;

# BinarySub
# to subtract two objects
export 
    BinarySub ::static := proc(self :: UnivariatePolynomialOverPowerSeriesObject,
                                other :: UnivariatePolynomialOverPowerSeriesObject,
                                $)
        local x := get_common_vname(self, other);
        local u := max(DEGREE(self), DEGREE(other));
        local l := min(DEGREE(self), DEGREE(other));
        local apoly := self:-upoly;
        local bpoly := other:-upoly;
        local new_upoly := Array(0 .. l, i -> PowerSeriesObject:-BinarySub(apoly[i], bpoly[i]));
        if DEGREE(self) < DEGREE(other) then
            new_upoly(u+1) := 0; # reallocate Array 
            for local i from l+1 to u do
                new_upoly[i] := PowerSeriesObject:-Negate(bpoly[i]);
            end do;
        elif DEGREE(self) > DEGREE(other) then
            new_upoly(u+1) := 0;
            for local i from l+1 to u do 
                new_upoly[i] := PowerSeriesObject:-DeepCopy(apoly[i]);
            end do;
        end if;
        return Object(UnivariatePolynomialOverPowerSeriesObject, new_upoly, x);
    end proc;

# Negate
# to negate the input object
export 
    Negate ::static := proc(self :: UnivariatePolynomialOverPowerSeriesObject, $)
        local new_upoly := map(PowerSeriesObject:-Negate, self:-upoly);
        return Object(UnivariatePolynomialOverPowerSeriesObject, new_upoly, self:-vname);
    end proc;

# BinaryMultiply
# to multiply two objects
export 
    BinaryMultiply ::static := proc(self :: UnivariatePolynomialOverPowerSeriesObject,
                                    other :: UnivariatePolynomialOverPowerSeriesObject,
                                    $)
        local x := get_common_vname(self, other);
        local d := DEGREE(self) + DEGREE(other);
        local new_upoly := Array(0 .. d);
        for local i from 0 to DEGREE(self) do
            for local j from 0 to DEGREE(other) do
                local product := PowerSeriesObject:-BinaryMultiply(self:-upoly[i], other:-upoly[j]);
                if new_upoly[i+j] = 0 then
                    new_upoly[i+j] := product;
                else 
                    new_upoly[i+j] := PowerSeriesObject:-BinaryAdd(new_upoly[i+j], product);
                end if;
            end do;
        end do;
        return Object(UnivariatePolynomialOverPowerSeriesObject, new_upoly, x);
    end proc;

# a local function used in NaryMultiply when size=constant 
# shouldn't be documented 
local
    multiply_binary_split :: static := proc(lterms :: nonemptylist(UnivariatePolynomialOverPowerSeriesObject), $)
        if numelems(lterms) = 1 then
            return lterms[1];
        else
            local k := trunc(numelems(lterms)/2);
            return BinaryMultiply(thisproc(lterms[1..k]), thisproc(lterms[k+1..]));
        end if;
    end proc;

# NaryMultiply
# to multiply a sequence of UPoPS objects 
# size : if the degree pattern is known, then use one of these three different modes 
#        {constant, increasing, decreasing} or leave it as default 
export
    NaryMultiply ::static := proc(factors :: seq({UnivariatePolynomialOverPowerSeriesObject, PowerSeriesObject, algebraic}),
                                  {size :: identical(constant, increasing, decreasing) := constant},
                                  $)
        local lfactors := [factors];
        if not type(lfactors, 'list'(UnivariatePolynomialOverPowerSeriesObject)) then
            local upops, rest;
            upops, rest := selectremove(type, lfactors, UnivariatePolynomialOverPowerSeriesObject);
            local vname := get_common_vname(op(upops));
            
            local pso := PowerSeriesObject:-NaryMultiply(op(rest));
            lfactors := [op(upops), pso:-ConvertToUPoPS(pso, vname)];
        end if;

        local n := numelems(lfactors);
        if n < 2 then
            return ifelse(n = 1, DeepCopy(lfactors[1]), One());
        elif size = constant then 
            # recursively spliting 
            return multiply_binary_split(lfactors);
        elif size = increasing then 
            # f(f(f(a, b), c), d)
            return foldl(BinaryMultiply, op(lfactors));
        else  
            # foldr(f, a, b, c, d) returns f(b, f(c, f(d, a))), so in order to have the smallest
            # entries nested the most deeply, we need to move the last entry to the front.
            return foldr(BinaryMultiply, lfactors[-1], op(.. -2, lfactors));
        end if;
    end proc;

export
    `*` ::static := proc(factors :: seq({UnivariatePolynomialOverPowerSeriesObject, PowerSeriesObject, algebraic}), $)
        return NaryMultiply(factors);
    end proc;

# Exponentiate
# to compute self^n 
# n : nonnegint 
export 
    Exponentiate ::static := proc(self :: UnivariatePolynomialOverPowerSeriesObject, 
                                n :: nonnegint, 
                                $)
        if n = 0 then 
            return One(); 
        elif n = 1 then 
            return DeepCopy(self);
        end if;        
        local p := self;
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

export
    `^` ::static := proc(x, y, $)
        if type(x, UnivariatePolynomialOverPowerSeriesObject) and type(y, ':-nonnegint') then
            return x:-Exponentiate(x, y);
        else 
            error "invalid input: this form of exponentiation is not implemented";
        end if;
    end proc;

# Note: Division algorithm for UPoPS is well-defined using 
# ConvertToPowerSeries and PSO:-Divide

#################
# ## Taylor Shift
#################

# generate an upper triangular matrix d+1 \times d+1 of (y + c)^i coefficients
# Note trying to use classic tricks to reduce cache misses in Maple causes slow down,
# to reduce memory usage, we ended up using only CoefficientList.
# used in taylorShift method 
local 
    taylor_shift_matrix ::static := proc(c :: {numeric, algebraic, algnum, algnumext, abstract_rootof, radical},
                                        d :: nonnegint,
                                        $)
        local mat := Matrix(d + 1, d + 1);
        local yc := 1; # will be (y+c)^(i-1)
        local isnum := type(c, 'numeric');
        local y;
        for local i from 1 to d+1 do
            local ycc := PolynomialTools:-CoefficientList(yc, y); # speed-up i > 90
            for local j from 1 to i do
                mat[j, i] := ycc[j];
            end do;
            yc := ifelse(isnum = true, expand(yc*(y + c)), Algebraic:-Expand(subsindets((yc*(y + c)), 'radical', convert, RootOf))); # AUTO_Expand is not efficient here! 
        end do;
        return mat;
    end proc;

# generate shifted upoly w.r.t 'mat', the upper triangular matrix from taylor_shift_matrix
# used in taylorShift method 
local 
    taylor_shift_make_upoly ::static := proc(self :: UnivariatePolynomialOverPowerSeriesObject,
                                                mat :: Matrix,
                                                $)
        local d := DEGREE(self);
        local new_upoly := Array(0 .. d);
        for local i from 0 to d do
            local j;
            new_upoly[i] := PowerSeriesObject:-NaryAdd(seq(self:-upoly[j], j = i .. d),
                                                       ':-coefficients' = [seq(mat[i+1, j+1], j = i .. d)]);
        end do;
        return new_upoly;
    end proc;

# Taylor Shift 
# f(y) -> f(y + c) where c is a numeric/algebraic.
export 
    TaylorShift ::static := proc(self :: UnivariatePolynomialOverPowerSeriesObject, 
                                c :: {numeric, algebraic, algnum, algnumext, abstract_rootof, radical},
                                $)
        if c = 0 then
            return self;
        end if;
        
        local mat := taylor_shift_matrix(c, DEGREE(self));
        return Object(UnivariatePolynomialOverPowerSeriesObject, taylor_shift_make_upoly(self, mat), self:-vname);
    end proc;
