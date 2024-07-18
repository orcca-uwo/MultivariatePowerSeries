export
Substitute :: static := proc(eqns :: {set, list, thistype}(name = {PowerSeriesObject, algebraic}),
                             _self :: PowerSeriesObject,
                             $)
local all_eqns := (if type(eqns, ':-list') then
                       eqns;
                   elif type(eqns, ':-set') then
                       [op(eqns)];
                   else
                       [eqns];
                   end if);
local vars := _self:-vars;
    all_eqns := select(eqn -> lhs(eqn) in vars, all_eqns);
    all_eqns := remove(evalb, all_eqns);
    if all_eqns = [] then
        return _self;
    elif numelems(convert(map(lhs, all_eqns), ':-set')) <> numelems(all_eqns) then
        error "these left hand sides occur more than once: %1", convert(ListTools:-FindRepetitions(
            map(lhs, all_eqns)), ':-set');
    end if;

local pso_eqns, poly_eqns, other_alg_eqns;
    pso_eqns, poly_eqns := selectremove(type, all_eqns, ':-anything' = PowerSeriesObject);
    poly_eqns, other_alg_eqns := selectremove(type, poly_eqns,
                                              ':-anything' = ':-polynom'(COEFFICIENT_TYPE));
local actually_polynomial;
    actually_polynomial, pso_eqns := selectremove(
        pso_eqn -> type(function:-GetAnalyticExpression(rhs(pso_eqn)),
                        ':-And'(':-polynom'(COEFFICIENT_TYPE), ':-Not'(undefined))), pso_eqns);
local eqn;
    poly_eqns := [op(poly_eqns),
                  seq(lhs(eqn) = function:-GetAnalyticExpression(rhs(eqn)),
                      eqn in actually_polynomial)];
local dummy;
    poly_eqns := map(eqn -> lhs(eqn) = AUTO_EXPAND(dummy, rhs(eqn)), poly_eqns);

    if pso_eqns = [] and type([op(poly_eqns), op(other_alg_eqns)], ':-nonemptylist'(
        ':-name' = '{:-`*`, :-thistype}'({':-name', ':-name'^':-rational'}))) then
        # Note that this excludes the constant case for the right hand side!
        poly_eqns := [op(poly_eqns), op(other_alg_eqns)];
        local x;
        return product_of_powers_case(_self, [op(poly_eqns), seq(x = x, x in vars minus convert(
            map(lhs, poly_eqns), ':-set'))]);

    elif other_alg_eqns <> [] then
        error "invalid substitution%{1|||s}, %2", numelems(other_alg_eqns), other_alg_eqns;
    end if;

local result := _self;
    
local shifts := table(_self:-vars =~ 0);
local any_nonzeroes := false;
    
    pso_eqns := [for local pso_eqn in pso_eqns do
                     local left := lhs(pso_eqn);
                     local right := rhs(pso_eqn);
#                         local (left, right) := op(pso_eqn);
                     local constant_coefficient := HOMOGENEOUS_PART(right, 0);
                     
                     if not Testzero(constant_coefficient) then
                         any_nonzeroes := true;
                         shifts[left] := constant_coefficient;
                         (* yield *) left = right - constant_coefficient;
                     else
                         (* yield *) pso_eqn;
                     end if;
                 end do];
    
    poly_eqns := [for local poly_eqn in poly_eqns do
                      local left := lhs(poly_eqn);
                      local right := rhs(poly_eqn);
#                          local (left, right) := op(poly_eqn);
                      local constant_coefficient := eval(
                          right, (indets(right, ':-name') minus {:-constants}) =~ 0);
                      
                      if not Testzero(constant_coefficient) then
                          any_nonzeroes := true;
                          shifts[left] := constant_coefficient;
                          (* yield *) left = right - constant_coefficient;
                      else
                          (* yield *) poly_eqn;
                      end if;
                  end do];
    
    if any_nonzeroes then
        try
            result := TaylorShift(result, [indices(shifts, ':-pairs')]);
        catch "invalid Taylor shift: tried to shift to a pole of the analytic expression":
            error "invalid substitution: tried to shift to a pole of the analytic expression";
        catch "cannot compute the Taylor shift of a power series of which the analytic expression is not known":
            error "substituting a unit power series into a power series with undefined analytic expression";
        end try;
    end if;
    
local zero_eqns;
    zero_eqns, poly_eqns := selectremove(type, poly_eqns, ':-name' = 0);
    
local renaming_eqns;
    renaming_eqns, poly_eqns := selectremove(type, poly_eqns, ':-name' = ':-name');

local poly_records := map(poly_eqn_to_record, poly_eqns);
local pso_records := map(pso_eqn_to_record, pso_eqns);

local all_records := [op(poly_records), op(pso_records)];
    all_records := sort(all_records, substitution_order);
local N := numelems(all_records);

    # If all_records[i]:-substituent contains x := all_records[j]:-substituee with i < j, then we
    # need to rename the occurrences of x in all_records[i]:-substituent to distinguish them from
    # the original occurrences x in _self, and then rename them back to x at the very end. Similar
    # if x is the left hand side of any renaming equation.
    for local i to N do
        local j;
        local to_rename := all_records[i]:-vars intersect {
            seq(all_records[j]:-substituee, j = i+1 .. N),
            seq(lhs(j), j in renaming_eqns)};
        all_records[i]:-renamings := map(v -> v = `tools/genglobal`(':-_X'), to_rename);
    end do;
local record;
local final_renamings := map(rhs = lhs, `union`(seq(record:-renamings, record in all_records))) union
    convert(renaming_eqns, ':-set');

local result_vars := result:-vars;

    if zero_eqns <> [] then
        result_vars := result_vars minus convert(map(lhs, zero_eqns), ':-set');
        result := Object(PowerSeriesObject, Array(0 .. 0, [result:-hpoly[0]]), 0,
                         subs_0_gen, result_vars,
                         ["parent" = result, "equations" = zero_eqns],
                         ifelse(type(result:-algexpr, undefined), undefined,
                                subs(zero_eqns, result:-algexpr)));
    end if;

local v;
    for record in all_records do
        result_vars := (result_vars minus {record:-substituee}) union subs(record:-renamings, record:-vars);

        if record:-type = "poly" then
            local ancestors := Record[record](
                "parent" = result,
                "terms" = convert_hpoly_from_poly(subs(record:-renamings, record:-substituent),
                                                  ':-output' = table),
                "degreelist",
                "termlist",
                "degree",
                "ldegree");
            local termlist := sort([indices(ancestors:-terms, ':-pairs')], ':-key' = lhs);
            ancestors:-degreelist := map(lhs, termlist);
            ancestors:-termlist := map(rhs, termlist);
            ancestors:-ldegree := ldegree(ancestors:-substituent);
            ancestors:-degree := degree(ancestors:-substituent);

            local result_algexpr := ifelse(type(result:-algexpr, undefined), undefined,
                                           eval(result:-algexpr, record:-substituee =
                                                subs(record:-renamings, record:-substituent)));

            result := Object(PowerSeriesObject, Array(0 .. 0, [result:-hpoly[0]]), 0,
                             subs_nonunit_poly_gen, result_vars,
                             ancestors,
                             result_algexpr);

        else
            ASSERT(record:-type = "pso");

            local ancestors := Record[record]("parent" = result);
            
            local result_algexpr := (if type(result:-algexpr, undefined)
                                     or type(record:-substituent:-algexpr, undefined) then
                                         undefined;
                                     else
                                         eval(result:-algexpr, record:-substituee =
                                              subs(record:-renamings, record:-substituent:-algexpr));
                                     end if);

            result := Object(PowerSeriesObject, Array(0 .. 0, [result:-hpoly[0]]), 0,
                             subs_nonunit_pso_gen, result_vars,
                             ancestors,
                             result_algexpr);
        end if;
    end do;

    if final_renamings <> {} then
        result_vars := subs(final_renamings, result_vars);
        result := Object(PowerSeriesObject, Array(0 .. 0, [result:-hpoly[0]]), 0,
                         renaming_gen, result_vars,
                         ["parent" = result, "eqns" = final_renamings],
                         subs(final_renamings, result:-algexpr));
    end if;

    return result;
end proc;

local
eqn_to_record :: static := (typ, hom0, vars) ->
proc(eqn :: equation, $)
local result := Record[':-packed'](
    "type"        = typ,
    "substituee"  = lhs(eqn),
    "substituent" = rhs(eqn),
    "renamings"   = {},
    "cache_table" = table(), # used for different things depending on the generator!
    "vars",
    "hom0");
    result:-hom0 := hom0(result:-substituent);
    result:-vars := vars(result:-substituent);

    return result;
end proc;

local
poly_eqn_to_record :: static := eqn_to_record("poly",
                                              poly -> eval(poly, indets(poly, ':-name') =~ 0),
                                              poly -> indets(poly, ':-name'));
local
pso_eqn_to_record :: static := eqn_to_record("pso",
                                             pso -> function:-HomogeneousPart(pso, 0),
                                             function:-Variables);

local
substitution_order :: static := proc(i, j)
    if i:-type = "pso" then
        if j:-type = "pso" then
            # Two psos. Keep them in their relative position.
            return true;
        else
            # Polys before psos.
            return false;
        end if;
    else
        if j:-type = "pso" then
            # Polys before psos.
            return true;
        else
            # Two polynomials. Sort by (degree - ldegree), then degree, then length.
            local i_properties := [degree, ldegree, length](i:-substituent);
            local j_properties := [degree, ldegree, length](j:-substituent);
            if i_properties[1] - i_properties[2] < j_properties[1] - j_properties[2] then
                return true;
            elif i_properties[1] - i_properties[2] > j_properties[1] - j_properties[2] then
                return false;
            elif i_properties[1] < j_properties[1] then
                return true;
            elif i_properties[1] > j_properties[1] then
                return false;
            elif i_properties[3] < j_properties[3] then
                return true;
            elif i_properties[3] > j_properties[3] then
                return false;
            end if;
            
            # Tied in everything. Keep current order.
            return true;
        end if;
    end if;
end proc;

local
renaming_gen :: static := proc(_self :: PowerSeriesObject, 
                               d :: nonnegint,
                               $)
    return subs(_self:-ancestors:-eqns, HOMOGENEOUS_PART(_self:-ancestors:-parent, d));
end proc;

local
subs_0_gen :: static := proc(_self :: PowerSeriesObject,
                             d :: posint,
                             $)
    return subs(_self:-ancestors:-equations,
                HOMOGENEOUS_PART(_self:-ancestors:-parent, d));
end proc;

# Let
#
#    f := %sum(a[i, n] * x^i * z^n, `i,n`=0..infinity)
#    g := %sum(b[j] * y^j, j=1..infinity)
#    fg := subs(x = g, f)
#
# then
#
#    fg = %sum(%sum(b[j] * y^j, j=1..infinity)^i * a[i, n] * z^n, `i,n`=0..infinity)
#
# {multinomial theorem, sum over (infinite) vectors alpha}
#
#       = %sum(%sum(multinomial(i, alpha) * prod((b[j] * y^j)^alpha[j], j=1..infinity),
#             %sum(alpha) = i) * a[i, n] * z^n, `i,n`=0..infinity)
#
# {substitute i = %sum(alpha); sum over *all* alpha (this is a swap of summation order}
#
#       = %sum(%sum(multinomial(%sum(alpha), alpha) * prod((b[j] * y^j)^alpha[j], j=1..infinity) *
#             a[%sum(alpha), n] * z^n, alpha), n=0..infinity)
#
# {group by %sum(j * alpha[j], j=1..infinity) =: m - n}
#
#       = %sum(%sum(%sum(multinomial(%sum(alpha), alpha) *
#               prod((b[j] * y^j)^alpha[j], j=1..infinity) * a[%sum(alpha), n] * z^n,
#               %sum(j * alpha[j], j=1..infinity) = m - n), m=n..infinity), n=0..infinity)
#
# {swap summation of m and n; no contributions with j > m-n}
#
#       = %sum(%sum(%sum(multinomial(%sum(alpha), alpha) *
#               prod((b[j] * y^j)^alpha[j], j=1..infinity) * a[%sum(alpha), n] * z^n,
#               %sum(j * alpha[j], j=1..m-n) = m - n), n=0..m), m=0..infinity)
#
# {select the homogeneous component of degree m}
#
#    fg[m] = %sum(%sum(multinomial(%sum(alpha), alpha) *
#                prod((b[j] * y^j)^alpha[j], j=1..infinity) * a[%sum(alpha), n] * z^n,
#                %sum(j * alpha[j], j=1..m-n) = m - n), n=0..m)
#
# {swap summation of n and alpha}
#
#          = %sum(%sum(multinomial(%sum(alpha), alpha) *
#                prod((b[j] * y^j)^alpha[j], j=1..infinity) * a[%sum(alpha), n] * z^n,
#                n=m-%sum(j*alpha[j], j=1..m)), %sum(j * alpha[j], j=1..m) = 0..m)
#
# {n has only a single value}
#
#          = %sum(multinomial(%sum(alpha), alpha) *
#              prod((b[j] * y^j)^alpha[j], j=1..infinity) *
#              a[%sum(alpha), m-%sum(j*alpha[j])] * z^(m-%sum(j*alpha[j])),
#              %sum(j * alpha[j], j=1..m) = 0..m).
#
# So: we find all alpha with that property, then collect all these terms.

# Possible future extension for multiple substitutions at once:
#
# Let
#
#   f := %sum(a[n] * x^n1 * y^n2 * z^n3, n = [0,0,0] .. [infinity,infinity,infinity])
#   g := %sum(b[j] * x^j, j=1..infinity)    (b[0] = 0)
#   h := %sum(c[j] * y^j, j=1..infinity)    (c[0] = 0)
#   fgh := subs({x = g, y = h}, f)
#
# then
#
#   fgh = %sum(a[n] * %sum(b[j] * x^j, j=1..infinity)^n1 * %sum(c[j] * y^j, j=1..infinity)^n2 *
#             z^n3, n = [0,0,0] .. [infinity,infinity,infinity])
#
# {multinomial theorem, sum over (infinite) vectors alpha, beta;
#  notation: |alpha| = %sum(alpha[j], j=1..infinity)}
#
#       = %sum(a[n] *
#             %sum(multinomial(|alpha|, alpha) * prod((b[j] * x^j)^alpha[j], j=1..infinity),
#                 |alpha| = n1) *
#             %sum(multinomial( |beta|,  beta) * prod((c[j] * y^j)^ beta[j], j=1..infinity),
#                 |beta|  = n2) * z^n3, n = [0,0,0] .. [infinity,infinity,infinity])
#
# {sum over *all* alpha, beta, instead of grouped by n2, n3 (this is a swap of summation order);
#  notation: |j alpha[j]| = %sum(j * alpha[j], j=1..infinity)}
#
#       = %sum(%sum(%sum(a[|alpha|, |beta|, n3] *
#             multinomial(|alpha|, alpha) * prod(b[j]^alpha[j], j=1..infinity) * x^|j alpha[j]| *
#             multinomial( |beta|,  beta) * prod(c[j]^ beta[j], j=1..infinity) * y^|j  beta[j]| *
#             z^n3, n3=0..infinity), alpha), beta)
#
# {introduce k = |j alpha[j]|, ell = |j beta[j]|, and group by k and j}
#
#       = %sum(%sum(%sum(%sum(%sum(a[|alpha|, |beta|, n3] *
#                 multinomial(|alpha|, alpha) * prod(b[j]^alpha[j], j=1..infinity) *
#                 multinomial( |beta|,  beta) * prod(c[j]^ beta[j], j=1..infinity) *
#                 x^k * y^ell * z^n3, |j alpha[j]| = k), |j beta[j]| = ell),
#             k = 0 .. infinity), ell = 0 .. infinity), n3 = 0 .. infinity)
#
# So: for a given homogeneous degree d, we find all alpha, beta, n3 such that |j alpha[j]| + |j
# beta[j]| + n3 = d and sum all corresponding terms. It should be possible to express that in a
# single weighted_vector_object... if we can make it keep track of |alpha| and |beta|
# separately. And it looks pretty clear how you would generalize this to arbitrary numbers of
# substitutions / non-substituted variables.

local
subs_nonunit_gen_worker :: static := proc(_self :: PowerSeriesObject,
                                          m :: posint,
                                          terms :: {table, Array},
                                          degreelist :: list(posint),
                                          $)
local ancestors := _self:-ancestors;
local result := Array([]);
local n_dl := numelems(degreelist);

local weighted_vector_object := Object(weighted_bounded_vectors, degreelist, m);
local alpha := weighted_vector_object:-Values();

    do
        local alphasum := weighted_vector_object:-Sum();
        # msjaj = m - %sum(j * alpha[j], j = 0 .. m).
        local msjaj := m - weighted_vector_object:-WeightSum();
    
        # We will be looking for terms a[alphasum,msjaj] * x^alphasum * z^msjaj, to be found in the
        # homogeneous component of degree alphasum + msjaj.
        local deg := alphasum + msjaj;

        # For these generators, cache_table caches the values of "by_x_degree".
        if not assigned(ancestors:-cache_table[deg]) then
            local hom_part := convert(HOMOGENEOUS_PART(ancestors:-parent, deg),
                                      ':-list', ':-`+`');
            # By substituting substituee=1, we turn a[alphasum,msjaj] * x^alphasum * z^msjaj
            # into a[alphasum,msjaj] * z^msjaj.
            # TODO: should this use PolynomialTools:-CoefficientVector?
            ancestors:-cache_table[deg] := table(
                ':-sparse',
                subs(ancestors:-substituee = 1,
                     map(eqn -> lhs(eqn) = add(rhs(eqn)), [entries(
                         ListTools:-Classify(term -> degree(term, ancestors:-substituee),
                                             hom_part), ':-pairs')])));
        end if;

        local by_x_degree := ancestors:-cache_table[deg];
        local zterms := by_x_degree[alphasum];
        if zterms <> 0 then
            local i;
            zterms *= mul(ifelse(alpha[i] = 0, 1, terms[degreelist[i]] ^ alpha[i]),
                          i = 1 .. n_dl);
            zterms *= combinat:-multinomial(alphasum,
                                            op(sort(subs(0 = NULL, [seq(alpha)]))));
            result ,= zterms;
        end if;
    until not weighted_vector_object:-Next();
    
local dummy;
    return AUTO_EXPAND(dummy, add(result));
end proc;

local
weighted_bounded_vectors :: static := module()
option object;
local
    state, out_alias;

local
    ModuleCopy :: static := proc(_self, other,
                                 weights :: list(posint),
                                 weightsum :: nonnegint,
                                 $)
    local n := numelems(weights);
        ASSERT(n = 0 or min(weights[2..] -~ weights[..-2]) >= 0);
        _self:-state := Vector([0,             # 1: current weight sum
                                0,             # 2: current sum
                                n,             # 3: n
                                weightsum,     # 4: requested max weightsum
                                op(weights),   # 5 .. n+4: weights
                                0 $ n,         # n+5 .. 2*n+4: values
                                NULL], ':-datatype' = ':-integer'[4]);
        _self:-out_alias := ArrayTools:-Alias(_self:-state, n+4, [n]);
    end proc;

$define CURRENT_WEIGHT_SUM state[1]
$define CURRENT_SUM        state[2]
$define N                  state[3]
$define REQ_WEIGHT_SUM     state[4]
$define WEIGHT(i)          state[4+i]
$define VALUE(i)           state[4+N+i]
    
export
    WeightSum :: static := proc(_self, $)
        return _self:-CURRENT_WEIGHT_SUM;
    end proc;

export
    Sum :: static := proc(_self, $)
        return _self:-CURRENT_SUM;
    end proc;

export
    Values :: static := proc(_self, $)
        return _self:-out_alias;
    end proc;

export
    Next :: static := proc(_self, $)
        return evalb(autocompiled_next(_self:-state) = 1);
    end proc;

local
    autocompiled_next :: static := proc(state :: Vector(datatype = integer[4]))
    option autocompile;
    local i :: integer := 1;
        while i <= N and CURRENT_WEIGHT_SUM + WEIGHT(i) > REQ_WEIGHT_SUM do
            # "Empty" values[i]
            CURRENT_WEIGHT_SUM := CURRENT_WEIGHT_SUM - VALUE(i) * WEIGHT(i);
            CURRENT_SUM := CURRENT_SUM - VALUE(i);
            VALUE(i) := 0;
            i := i + 1;
        end do;

        if i > N then
            return 0;
        end if;

        # Increase values[i] by 1
        CURRENT_WEIGHT_SUM := CURRENT_WEIGHT_SUM + WEIGHT(i);
        CURRENT_SUM := CURRENT_SUM + 1;
        VALUE(i) := VALUE(i) + 1;

        return 1;
    end proc;

$undef CURRENT_WEIGHT_SUM
$undef CURRENT_SUM
$undef N
$undef REQ_WEIGHT_SUM
$undef WEIGHT
$undef VALUE
end module;

local
subs_nonunit_poly_gen :: static := proc(_self :: PowerSeriesObject, 
                                        d :: posint,
                                        $)
    return subs_nonunit_gen_worker(_self, d, _self:-ancestors:-terms,
                                   select(`<=`, _self:-ancestors:-degreelist, d));
end proc;

local
subs_nonunit_pso_gen :: static := proc(_self :: PowerSeriesObject, 
                                       d :: nonnegint,
                                       $)
local substituent := _self:-ancestors:-substituent;
    substituent:-ensure_degree(substituent, d);

    # nonunit, so we know there's no degree 0 term; and we need to do this with seq so that we can
    # access hpoly directly. Note that the renamings still need to be applied, whereas for the poly
    # case, we can do this beforehand.
local i;
local degreelist := [seq(ifelse(substituent:-hpoly[i] = 0, NULL, i), i = 1 .. d)];
local terms := table(subs(_self:-ancestors:-renamings, [seq(i = substituent:-hpoly[i], i in degreelist)]));

    return subs_nonunit_gen_worker(_self, d, terms, degreelist);
end proc;

(*

This part is obsolete with the introduction of Taylor shifts - feel free to remove

# Let
#
#    f := %sum(a[i, n] * x^i * z^n, `i,n`=0..infinity)
#    g := %sum(b[j] * y^j, j=0..infinity)
#    fg := subs(x = g, f)
#
# then (see MPL-23990)
#
#    fg = %sum(a[i, n] * %sum(b[j] * y^j, j=0..infinity)^i * z^n, `i,n`=0..infinity)
#
# = { split off b[0] }
#
#         %sum(a[i, n] * (b[0] + %sum(b[j] * y^j, j=1..infinity))^i * z^n, `i,n`=0..infinity)
#
# = { binomial theorem }
#
#         %sum(a[i, n] * %sum(binomial(i, k) * b[0]^(i-k) * %sum(b[j] * y^j, j=1..infinity)^k,
#             k=0..i) * z^n, `i,n`=0..infinity)
#
# = { multinomial theorem, sum over (infinite) vectors alpha }
#
#         %sum(a[i, n] * %sum(binomial(i, k) * b[0]^(i-k) *
#             %sum(multinomial(k, alpha) * prod((b[j] * y^j)^alpha[j], j=1..infinity),
#               %sum(alpha) = k), k=0..i) * z^n, `i,n`=0..infinity)
#
# = { bring all factors into sum over alpha }
#
#         %sum(%sum(%sum(a[i, n] * binomial(i, k) * b[0]^(i-k) *
#               multinomial(k, alpha) * prod((b[j] * y^j)^alpha[j], j=1..infinity) * z^n,
#               %sum(alpha) = k), k=0..i), `i,n`=0..infinity)
#
# = { swap summation over k and i (changing their ranges), then swap i and alpha }
#
#         %sum(%sum(%sum(a[i, n] * binomial(i, k) * b[0]^(i-k) *
#               multinomial(k, alpha) * prod((b[j] * y^j)^alpha[j], j=1..infinity) * z^n,
#               i=k..infinity), %sum(alpha) = k), `k,n`=0..infinity)
#
# = { take factors independent of i out of inner sum }
#
#         %sum(%sum(%sum(a[i, n] * binomial(i, k) * b[0]^(i-k), i=k..infinity) *
#             multinomial(k, alpha) * prod((b[j] * y^j)^alpha[j], j=1..infinity) * z^n,
#             %sum(alpha) = k), `k,n`=0..infinity)
#
# = { the sums over k, n, and alpha iterate over all (nonnegative integer n, sparse nonnegative
#     integer vector alpha), grouped by k = %sum(alpha); we can write that instead as a sum over all
#     (nonnegative integer n, nonnegative integer vector alpha), grouped by m = n + %sum(j*alpha[j],
#     j=1..infinity); the sum with range %sum(j*alpha[j])=m-n is over alpha, in case it's unclear }
#
#         %sum(%sum(%sum(%sum(a[i, n] * binomial(i, %sum(alpha)) * b[0]^(i-%sum(alpha)),
#                 i=%sum(alpha)..infinity) *
#               multinomial(%sum(alpha), alpha) *
#               %prod((b[j] * y^j)^alpha[j], j=1..infinity) * z^n,
#               %sum(j*alpha[j])=m-n), n=0..m), m=0..infinity)
#
# The homogeneous degree of each term in the sum over m, is m, so we have:
#
#    fg[m] = %sum(%sum(%sum(a[i, n] * binomial(i, %sum(alpha)) * b[0]^(i-%sum(alpha)),
#                  i=%sum(alpha)..infinity) *
#                multinomial(%sum(alpha), alpha) *
#                %prod((b[j] * y^j)^alpha[j], j=1..infinity) * z^n,
#                %sum(j*alpha[j])=m-n), n=0..m)
#
# This expression for the m'th homogeneous component still contains infinite sums. To deal with
# these, we use derivatives.
#
#         %diff(f, x $ k)
#
# = { definition of f, polynomial derivatives }
#
#         %sum(%sum(a[i, n] * i!/(i-k)! * x^(i-k), i=k..infinity) * z^n, n=0..infinity)
#
# Dividing by k! and substituting k = %sum(alpha), x = b[0], we get:
#
#         %eval(%diff(f, x $ %sum(alpha)) / (%sum(alpha))!, x = b[0])
#
# = { definition of binomial }
#
#         %sum(%sum(a[i, n] * binomial(i, %sum(alpha)) * b[0]^(i-%sum(alpha)),
#             i=%sum(alpha)..infinity) * z^n, n=0..infinity)
#
# We take the derivative w.r.t. z, h times.
#
#         %eval(%diff(f, x $ %sum(alpha), z $ h) / (%sum(alpha))!, x = b[0])
#
# = { definition of derivative }
#
#         %sum(%sum(a[i, n] * binomial(i, %sum(alpha)) * b[0]^(i-%sum(alpha)),
#             i=%sum(alpha)..infinity) * %diff(z^n, z $ h), n=0..infinity)
#
# Evaluate at z = 0.
#
#         %eval(%diff(f, x $ %sum(alpha), z $ h) / (%sum(alpha))!, {x = b[0], z = 0})
#
# = { the only term with a nonzero contribution has n=h }
#
#         %sum(a[i, h] * binomial(i, %sum(alpha)) * b[0]^(i-%sum(alpha)),
#           i=%sum(alpha)..infinity) * h!
#
# Rename h=n and divide by n!.
#
#         %eval(%diff(f, x $ %sum(alpha), z $ n) / (%sum(alpha))! / n!, {x = b[0], z = 0})
#
# = { elementary algebra }
#
#         %sum(a[i, n] * binomial(i, %sum(alpha)) * b[0]^(i-%sum(alpha)),
#           i=%sum(alpha)..infinity)
#
# Substitute this into the expression for fg[m] above.
#
#    fg[m] = %sum(%sum(%eval(%diff(f, x $ %sum(alpha), z $ n) / (%sum(alpha))! / n!,
#                  {x = b[0], z = 0}) *
#                multinomial(%sum(alpha), alpha) *
#                %prod((b[j] * y^j)^alpha[j], j=1..infinity) * z^n,
#                %sum(j*alpha[j])=m-n), n=0..m)
#
# What does this look like if there are multiple z's? I.e., redefine
#
#    f := %sum(%sum(a[i, n, l] * x^i * z^(n - l) * w^l, l=0..n), `i,n`=0..infinity)
#
# (with n indicating the homogeneous degree *in w and z*, instead of just degree in z). Then
#
#    fg[m] = %sum(%sum(%sum(%sum(a[i, n, l] * binomial(i, %sum(alpha)) * b[0]^(i-%sum(alpha)),
#                    i=%sum(alpha)..infinity) *
#                  multinomial(%sum(alpha), alpha) *
#                  %prod((b[j] * y^j)^alpha[j], j=1..infinity) * z^(n-l) * w^l,
#                  l=0..n), %sum(j*alpha[j])=m-n), n=0..m).
#
# We would still have
#
#         %eval(%diff(f, x $ %sum(alpha)) / (%sum(alpha))!, x = b[0])
#
# = { see above }
#
#         %sum(%sum(%sum(a[i, n, l] * binomial(i, %sum(alpha)) * b[0]^(i-%sum(alpha)),
#             i=%sum(alpha)..infinity) * z^(n-l) * w^l, l=0..n), n=0..infinity)
#
# but now we would differentiate w.r.t. z, hz times, and w.r.t. w, hw times:
#
#         %eval(%diff(f, x $ %sum(alpha), z $ hz, w $ hw) / (%sum(alpha))!, x = b[0])
#
# = { differentiation commutes with summation and constant multiplication }
#
#         %sum(%sum(%sum(a[i, n, l] * binomial(i, %sum(alpha)) * b[0]^(i-%sum(alpha)),
#             i=%sum(alpha)..infinity) * diff(z^(n-l), z $ hz) * diff(w^l, w $ hw), l=0..n),
#           n=0..infinity)
#
# We evaluate at {w=0, z=0} (and assume 0 <= hz <= n-l, 0 <= hw <= l).
#
#         %eval(%diff(f, x $ %sum(alpha), z $ hz, w $ hw) / (%sum(alpha))!,
#           {x = b[0], w = 0, z = 0})
#
# = { only one term contributes to two of the sums: l = hw, n = hz + hw (the latter so that n-l = hz) }
#
#         %sum(a[i, hz + hw, hw] * binomial(i, %sum(alpha)) * b[0]^(i-%sum(alpha)),
#             i=%sum(alpha)..infinity) * hz! * hw!
#
# We rename hw = l, hz = n-l, and divide by l! * (n-l)!.
#
#         %eval(%diff(f, x $ %sum(alpha), z $ (n-l), w $ l) / (%sum(alpha))!,
#           {x = b[0], w = 0, z = 0}) / l! / (n-l)!
#
# = { algebra }
#
#         %sum(a[i, n, l] * binomial(i, %sum(alpha)) * b[0]^(i-%sum(alpha)),
#             i=%sum(alpha)..infinity)
#
# and now we're at a point again where we can substitute this into fg[m]:
#
#    fg[m] = %sum(%sum(%sum(%eval(%diff(f, x $ %sum(alpha), z $ (n-l), w $ l) / (%sum(alpha))!,
#                    {x = b[0], w = 0, z = 0}) / l! / (n - l)! *
#                  multinomial(%sum(alpha), alpha) *
#                  %prod((b[j] * y^j)^alpha[j], j=1..infinity) * z^(n-l) * w^l,
#                  l=0..n), %sum(j*alpha[j])=m-n), n=0..m)
#
# = { the outer two sums iterate over all alpha with a weighted sum of at most m, and have n = m -
#     %sum(j*alpha[j]); rewrite it as such explicitly and substitute for n }
#
#            %sum(%sum(%eval(%diff(f, x $ %sum(alpha), z $ (m-%sum(j*alpha[j])-l), w $ l) / (%sum(alpha))!,
#                  {x = b[0], w = 0, z = 0}) / l! / (m - %sum(j*alpha[j]) - l)! *
#                multinomial(%sum(alpha), alpha) *
#                %prod((b[j] * y^j)^alpha[j], j=1..infinity) * z^(m - %sum(j*alpha[j])-l) * w^l,
#                l=0..m-%sum(j*alpha[j])), %sum(j*alpha[j]) <= m)
#
# = { expand the multinomial and cancel against %sum(alpha)! }
#
#            %sum(%sum(%eval(%diff(f, x $ %sum(alpha), z $ (m-%sum(j*alpha[j])-l), w $ l),
#                  {x = b[0], w = 0, z = 0}) / l! / (m - %sum(j*alpha[j]) - l)! *
#                (1 / %prod(alpha[j]!)) *
#                %prod((b[j] * y^j)^alpha[j], j=1..infinity) * z^(m - %sum(j*alpha[j])-l) * w^l,
#                l=0..m-%sum(j*alpha[j])), %sum(j*alpha[j]) <= m)
#
# This indicates how the general case works: for each alpha, we differentiate f with respect to x $
# %sum(alpha), and with respect to all other variables v some number of times t[v] such that the sum
# of the t-values is m-%sum(j*alpha[j]). We find the product of (b[j] * y^j)^alpha[j], where if g is
# multivariate, b[j] * y^j represents the homogeneous component of degree j; multiply this by the
# appropriate constant factors; and we're done. Then sum over all such t-assignments and over all
# alpha.

local
subs_unit_gen_worker :: static := proc(_self :: PowerSeriesObject,
                                       m :: posint,
                                       terms :: {table, Array},
                                       degreelist :: list(posint),
                                       $)
local ancestors := _self:-ancestors;
local result := Array([]);
local n_dl := numelems(degreelist);

local vars := convert(ancestors:-parent:-vars minus {ancestors:-substituee}, ':-list');
local nvars := numelems(vars);
    
local weighted_vector_object := Object(weighted_bounded_vectors, degreelist, m);
local alpha := weighted_vector_object:-Values();

    do
        local alphasum := weighted_vector_object:-Sum();
        # msjaj = m - %sum(j * alpha[j], j = 0 .. m).
        local msjaj := m - weighted_vector_object:-WeightSum();

        # bjyjaj is the product of the factors (b[j] * y^j)^alpha[j], with corresponding constant
        # factors
        local i;
        local bjyjaj := mul(ifelse(alpha[i] = 0, 1, terms[degreelist[i]] ^ alpha[i]),
                            i = 1 .. n_dl);
        bjyjaj /= map[':-fold' = (:-`*`, 1)](factorial, alpha);

        for local tvalues in ifelse(
            nvars = 1, [[msjaj]], Iterator:-BoundedComposition([msjaj $ nvars], msjaj)) do

            # For these generators, cache_table caches the values of the derivatives, multiplied
            # appropriately.
            local tbl_index := [alphasum, seq(tvalues)];
            if not assigned(ancestors:-cache_table[op(tbl_index)]) then
                if assign_unit_gen_worker_cache_entry(ancestors, tbl_index, vars) = 0 then
                    next;
                end if;
                # maybe if any derivative is identically zero, we can somehow cut something off in
                # tvalues?
            end if;

            result ,= bjyjaj * ancestors:-cache_table[op(tbl_index)];
        end do;
    until not weighted_vector_object:-Next();

local dummy;
    return AUTO_EXPAND(dummy, add(result));
end proc;

local
assign_unit_gen_worker_cache_entry :: static := proc(ancestors :: record,
                                                     tbl_index :: list,
                                                     vars :: list,
                                                     $)
local i;
local deriv := diff(ancestors:-parent:-algexpr, ancestors:-substituee $ tbl_index[1],
                    seq(vars[i] $ tbl_index[i+1], i = 1 .. nops(vars)));
    deriv := eval(deriv, {ancestors:-substituee = ancestors:-hom0, seq(vars =~ 0)});
    if deriv = 0 then
        ancestors:-cache_table[op(tbl_index)] := 0;
    else
        ancestors:-cache_table[op(tbl_index)] := deriv * mul(vars[i] ^ tbl_index[i+1] / tbl_index[i+1]!, i = 1 .. nops(vars));
    end if;
    return ancestors:-cache_table[op(tbl_index)];
end proc;

local
subs_unit_poly_gen :: static := proc(_self :: PowerSeriesObject, 
                                     d :: nonnegint,
                                     $)
    return subs_unit_gen_worker(_self, d, _self:-ancestors:-terms,
                                subs(0 = NULL, select(`<=`, _self:-ancestors:-degreelist, d)));
end proc;

local
subs_unit_pso_gen :: static := proc(_self :: PowerSeriesObject, 
                                    d :: nonnegint,
                                    $)
local substituent := _self:-ancestors:-substituent;
    substituent:-ensure_degree(substituent, d);

    # the degree 0 term is taken care of separately; and we need to do this with seq so that we can
    # access hpoly directly. Note that the renamings still need to be applied, whereas for the poly
    # case, we can do this beforehand.
local i;
local degreelist := [seq(ifelse(substituent:-hpoly[i] = 0, NULL, i), i = 1 .. d)];
local terms := table(subs(_self:-ancestors:-renamings, [seq(i = substituent:-hpoly[i], i in degreelist)]));

    return subs_unit_gen_worker(_self, d, terms, degreelist);
end proc;
*)

local
product_of_powers_case :: static := proc(_self :: PowerSeriesObject,
                                         eqns :: list(name = algebraic),
                                         $)
local vars_in := map(lhs, eqns);
local n := numelems(vars_in);
local vars_out := convert(indets(map(rhs, eqns), ':-name'), ':-list');
local vars_out_table := table([for local i, nm in vars_out do nm = i end do]);
local m := numelems(vars_out);

local as_products := map(eqn -> convert(rhs(eqn), ':-list', ':-`*`'), eqns);
local as_powers := map2(map, convert, as_products, ':-list', ':-`^`');
local homogeneous_degrees := map(lst -> add(map(`?[]`, lst, [2])), as_powers);

local mindeg := min(homogeneous_degrees);
    if mindeg <= 0 then
        local i := min[':-index'](homogeneous_degrees);
        error "unexpected substitution with nonpositive total degree: %1", eqns[i];
    end if;
local maxdeg := max(homogeneous_degrees);

    return Object(PowerSeriesObject, Array(0 .. 0, [_self:-hpoly[0]]), 0,
                  product_of_powers_subs_gen, convert(vars_out, ':-set'),
                  ["parent" = _self,
                   "term_cache" = table(), "mindeg" = mindeg, "maxdeg" = maxdeg,
                   "homogeneous_degrees" = homogeneous_degrees,
                   "substitutions" = eqns,
                   "vars_in" = vars_in, "vars_out" = vars_out,
                   NULL],
                  subs(eqns, _self:-algexpr));
end proc;

local
product_of_powers_subs_gen :: static := proc(_self :: PowerSeriesObject,
                                             d :: nonnegint,
                                             $)
local ancestors := _self:-ancestors;
local result := 0;

    for local deg from ceil(d / ancestors:-maxdeg) to floor(d / ancestors:-mindeg) do
        if not assigned(ancestors:-term_cache[deg]) then
            local all_terms := subs(ancestors:-substitutions, function:-HomogeneousPart(
                ancestors:-parent, deg));
            local as_sum := convert(all_terms, ':-list', ':-`+`');
            ancestors:-term_cache[deg] := ListTools:-Classify(degree, as_sum);
            if assigned(ancestors:-term_cache[deg][FAIL]) then
                error "non-polynomial result for power series: %1",
                ancestors:-term_cache[deg][FAIL];
            end if;
        end if;

        local terms := ancestors:-term_cache[deg][d];
        if type(terms, ':-set') then
            result += add(terms);
        end if;
    end do;

    if not type(result, ':-polynom') then
        error "non-polynomial result for power series: %1", result;
    end if;
    
    return result;
end proc;
