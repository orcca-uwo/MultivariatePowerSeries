###################################
## # Weierstrass Preparation UPoPS
###################################

# weierstrassUpdate
# given p = y^d + \sum_{i=0}^{d-1} b_i y^i 
# given alpha = \sum_{i=0}^{m} c_i y^i
# given F s.t. F = Array(F_i, i = 0, ..., d-1) where F_i = \sum_{k=0}^{i-1} b_k c_{i-k} 
# given r s.t. b_0, ..., b_{d-1}, c_0, ..., c_m known modulo M^r
# return b_0, ..., b_{d-1}, c_0, ..., c_m known module M^{r+1} in-place.
local 
    weierstrassUpdate ::static := proc(p :: UnivariatePolynomialOverPowerSeriesObject,
                                alpha :: UnivariatePolynomialOverPowerSeriesObject,
                                F :: Array,
                                r :: nonnegint,
                                $)
        local d := DEGREE(p);
        local s, dummy, k;
            # update b_0, ..., b_{d-1} module M^{r+1}
            for local i from 0 to d-1 do 
                s := add(seq(AUTO_EXPAND(dummy, PowerSeriesObject:-HomogeneousPart_local(p:-upoly[i], r-k) * PowerSeriesObject:-HomogeneousPart_local(alpha:-upoly[0], k)), k = 1 .. r-1));
                s := PowerSeriesObject:-HomogeneousPart_local(F[i], r) - s;
                PowerSeriesObject:-SetHomogeneousPart(p:-upoly[i], r, AUTO_EXPAND(dummy, s/HomogeneousPart_local(alpha:-upoly[0], 0)));
            end do;
    end proc;

# weierstrassGenerator
# given p = y^d + \sum_{i=0}^{d-1} b_i y^i 
# given alpha = \sum_{i=0}^{m} c_i y^i
# given F s.t. F = Array(F_i, i = 0, ..., d-1) where F_i = \sum_{k=0}^{i-1} b_k c_{i-k} 
# given requestedDeg : the requested deg from wp_gen, and requestedPSDeg : the index of the PowerSeries in upoly.
# update all PowerSeries (coefficients of) p and alpha up to degree requestedDeg
# return the updated PowerSeries at index requestedPSDeg in p. 
local 
    weierstrassGenerator ::static := proc(p :: UnivariatePolynomialOverPowerSeriesObject,
                                    alpha :: UnivariatePolynomialOverPowerSeriesObject,
                                    F :: Array,
                                    requestedDeg :: nonnegint,
                                    requestedPSDeg :: nonnegint := undefined,
                                    $)
        # print("In wGenerator, upperbound", upperbound(p:-upoly[0]:-hpoly));
        local d := DEGREE(p);
        if requestedPSDeg <> undefined and requestedPSDeg > d then
            error "requestedPSDeg is out of index.";
        elif UnivariatePolynomialOverPowerSeriesObject:-ApproximatelyZero(p, -2) then
            return p;
        end if;
        local u;
        local minDeg := min(seq(PowerSeriesObject:-GetPrecision(u), u in p:-upoly)); 
        for local r from minDeg+1 to requestedDeg do
            weierstrassUpdate(p, alpha, F, r);
        end do;
        if requestedPSDeg <> undefined then 
            return PowerSeriesObject:-HomogeneousPart_local(p:-upoly[requestedPSDeg], requestedDeg);
        end if;
    end proc;

# pDegree
# check if the input UPoPS is a valid input in weierstrassPreparation
# return the degree of p in weierstrassPreparation
local 
    pDegree ::static := proc(up :: UnivariatePolynomialOverPowerSeriesObject, $)
        if UnivariatePolynomialOverPowerSeriesObject:-ApproximatelyZero(up, 0) then
            error "invalid input: expected univariate polynomial over power series with at least "
            "one coefficient that is a unit, but received one that has no unit coefficients";
        end if;

        for local i from 0 to DEGREE(up) do
            if PowerSeriesObject:-IsUnit(up:-upoly[i]) then
                return i; # break;
            end if;
        end do;

        error "invalid input: expected univariate polynomial over power series with at least "
        "one coefficient that is a unit, but received one that has no unit coefficients";
    end proc;

# wp_gen
# wp generator 
local 
    wp_gen ::static := proc(self :: PowerSeriesObject, 
                    d :: nonnegint,
                    $)
        return weierstrassGenerator(PowerSeriesObject:-GetAncestor(self, 'p'), PowerSeriesObject:-GetAncestor(self, 'alpha'),
                                    PowerSeriesObject:-GetAncestor(self, 'F'), d, PowerSeriesObject:-GetAncestor(self, 'PSDeg'));
    end proc;

# WeierstrassPreparation
# given up an arbitrary UPoPS
# reutrn a unique pair (p, alpha) s.t. 
# p is a monic polynomial of degree d (p = y^d + \sum_{i=0}^{d-1} b_i y^i),
# alpha is an invertible power series (alpha = \sum_{i=0}^{m} c_i y^i),
# and up = alpha * p holds
# proc : if defined, then return the results of UpdatePrecision(alpha, proc) and UpdatePrecision(p, prec) in a optimized way. 
export 
    WeierstrassPreparation ::static := proc(up :: UnivariatePolynomialOverPowerSeriesObject, 
                                            prec :: nonnegint := undefined, 
                                            $)
        local d := pDegree(up);
        if d = 0 then
            return UnivariatePolynomialOverPowerSeriesObject:-One(), DeepCopy(up);
        end if;

        local m := DEGREE(up) - d;
        local F_out := Array(0 .. d-1); 
        local p_out, alpha_out;
        if up:-vname <> undefined then 
            p_out := Object(UnivariatePolynomialOverPowerSeriesObject, Array(0 .. d), up:-vname, ':-checkvariables' = false);
            alpha_out := Object(UnivariatePolynomialOverPowerSeriesObject, Array(0 .. m), up:-vname, ':-checkvariables' = false);
        else 
            p_out := Object(UnivariatePolynomialOverPowerSeriesObject, Array(0 .. d), ':-checkvariables' = false);
            alpha_out := Object(UnivariatePolynomialOverPowerSeriesObject, Array(0 .. m), ':-checkvariables' = false);
        end if;
        local b := p_out:-upoly; 
        local c := alpha_out:-upoly;

        local pso, all_variables := `union`(seq(pso:-Variables(pso), pso in up:-upoly));
        # initializing b 
        for local i from 0 to d-1 do 
            b[i] := Object(PowerSeriesObject, Array(0 .. 0, [PowerSeriesObject:-HomogeneousPart(up:-upoly[i], 0)]), 0,
                           wp_gen, all_variables, ['p' = p_out, 'alpha' = alpha_out, 'F' = F_out, 'PSDeg' = i]);
        end do;
        b[d] := PowerSeriesObject:-One();

        # initializing c
        c[m] := up:-upoly[d+m];
        for local j from m-1 to 0 by -1 do 
            c[j] := up:-upoly[d+j];
            for local i from d-1 to 0 by -1 do
                local k := d + j - i;
                if k <= m then
                    c[j] := PowerSeriesObject:-BinarySub(c[j], PowerSeriesObject:-BinaryMultiply(b[i], c[k]));
                end if;
            end do;
        end do;

        # initializing F_out
        F_out[0] := up:-upoly[0];
        for local j from 1 to d-1 do
            F_out[j] := up:-upoly[j];
            for local i from 1 to min(j, m) do
                F_out[j] := PowerSeriesObject:-BinarySub(F_out[j], PowerSeriesObject:-BinaryMultiply(b[j - i], c[i]));
            end do;
        end do;

        if prec <> undefined then
            weierstrassGenerator(PowerSeriesObject:-GetAncestor(p_out:-upoly[0], 'p'), 
                                 PowerSeriesObject:-GetAncestor(p_out:-upoly[0], 'alpha'), 
                                 PowerSeriesObject:-GetAncestor(p_out:-upoly[0], 'F'), 
                                 prec);
        end if;
        return p_out, alpha_out;
    end proc;
