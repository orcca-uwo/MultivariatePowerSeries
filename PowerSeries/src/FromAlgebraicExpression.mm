export
FromAlgebraicExpression :: static := proc(expr :: algebraic,
                                          {method :: identical(diff, series) := series},
                                          {force :: truefalse := false},
                                          {expand :: truefalse := true},
                                          $)
local variables := indets(expr, ':-name') minus {:-constants};
    variables := select[2](depends, expr, variables);
    if (not force) and type(expr, ':-ratpoly'(COEFFICIENT_TYPE, variables)) then
        return PowerSeries(expr);
    end if;

    if numelems(variables) > 1 and expand then
        local result := try_expand(expr, variables, _options['method']);
        if result <> FAIL then
            return result;
        end if;
    end if;

local hom0 := (try
                   Array(0 .. 0, [eval(expr, variables =~ 0)]);
               catch "numeric exception: division by zero":
                   error "could not form power series: the expression %1 appears to have a pole at zero", expr;
               end try);
   
    if method = diff then
        local derivatives := table([(0 $ numelems(variables)) = expr]);
        
        return Object(PowerSeriesObject,
                      hom0,
                      0,
                      algebraic_expression_gen_diff,
                      variables,
                      ["derivatives" = derivatives,
                       "variables" = convert(variables, ':-list'),
                       "max_derivatives_computed" = 0,
                       NULL],
                      expr,
                      NULL);
    else
        local used_order := table(variables =~ Order);
        
        local result := Object(PowerSeriesObject,
                               hom0,
                               0,
                               algebraic_expression_gen_series,
                               variables,
                               ["used_order" = used_order,
                                "found_order" = table(),
                                "degrees_stored" = 0,
                                NULL],
                               expr,
                               NULL);
        do_update_series(result);
        return result;
    end if;
end proc;

local
try_expand :: static := module()
local
    ModuleApply := proc(expr,
                        variables,
                        {method :: identical(diff, series) := series},
                        $)
        # If we can split the expression into multiple subexpressions, combined with algebraic
        # operations, that'll be much more efficient. E.g. exp(x+y+z) -> exp(x) * exp(y) * exp(z) or
        # even sin(x*y + y*z) -> sin(z*y)*cos(x*y)+cos(z*y)*sin(x*y).
        #
        # TODO: we would like to keep sin(x + x^2 + y + y^2) expanded only partly: as
        # sin(x^2+x)*cos(y+y^2) + cos(x^2+x)*sin(y+y^2). I'm not sure how to do that.
        #
        # TODO: avoid expanding (x + y)^100000 or sin(add(x[i], i=1..100)).
    local my_expr := expand(expr);
    local top_level_non_algebraic := indets[':-flat'](my_expr, ':-Not'(
        {':-`+`', ':-`*`', ':-`constant`', ':-anything'^':-posint'}));

    local variable_sets := map(na -> indets(na, ':-name') minus {:-constants}, top_level_non_algebraic);
    local max_variables := max(map(numelems, variable_sets));
        if max_variables >= numelems(variables) then
            return FAIL;
        end if;
        
        # Now we want to group children of the same DAG that have the same set of variables, or
        # where their set of variables is a subset of the others; e.g., for
        #
        #   sin(x) * sin(x^2*y) * sin(y) * z,
        #
        # we'd like to use
        #
        #   FromAlgebraicExpression(sin(x) * sin(x^2*y) * sin(y)) * FromAlgebraicExpression(z).
        #
        # On the other hand, we'd like to use one power series object if the same subexpression
        # occurs twice, in, say, y*sin(x) + z*sin(x). The following is not perfect, but a reasonable
        # approximation: we mark maximal subexpressions with a function call as in the first
        # example, and then create substitutions for them that should catch any duplicates.
    local from_algebraic_subexpression;
    local result := build_from_subexpressions(my_expr, from_algebraic_subexpression);

        # We should really use evalindets here, but MPL-26408 prevents it.
    local marked_subexpressions := indets(result, ':-specfunc'(from_algebraic_subexpression));
        return eval(result, map(expr -> expr = FromAlgebraicExpression(op(expr), _rest),
                                marked_subexpressions, _options['method'], ':-expand' = false));
    end proc;

local
    build_from_subexpressions := proc(expr,
                                      marker,
                                      $)
        if type(expr, {':-`*`', ':-`+`'}) then
            local i, t;
            local N := nops(expr);
            # A list, the ith operand of which is the finest set of sets of names into which op(i,
            # expr) can be split.
            local finest_variable_sets_per_op := [seq(get_finest_variable_sets(t), t in expr)];
            # A set of sets of names that are the finest sets of names that each subexpression can
            # be split into.
            local finest_variable_sets := `union`(op(finest_variable_sets_per_op));
            if `union`(op(finest_variable_sets)) in finest_variable_sets then
                # All variables necessarily occur together somewhere, so we should create the
                # expression as a single object.
                return marker(expr);
            end if;

            # i in splittable iff it might make sense to split op(i, expr) into different objects.
            local (nonsplittable, splittable) := selectremove(
                i -> member(`union`(op(finest_variable_sets_per_op[i])),
                            finest_variable_sets_per_op[i]), [seq(1 .. N)]);
            local nonsplittable_sets := {seq(op(finest_variable_sets_per_op[i]),
                                             i in nonsplittable)};
            local maximal_nonsplittable := remove(
                s -> ormap(t -> s subset t and s <> t, nonsplittable_sets), nonsplittable_sets);

            # Any operand whose set of names is contained in a maximal_nonsplittable can be grouped
            # together with other such operands "for free".
            local variable_sets_per_op := map(s -> `union`(op(s)), finest_variable_sets_per_op);
            local grouped := ListTools:-Classify(classifier, [seq(1 .. N)],
                                                 variable_sets_per_op, maximal_nonsplittable);

            # For any operands whose set of names is *not* contained in a maximal_nonsplittable,
            # i.e., op(i, expr) where i in grouped["none"], we recurse. The others are grouped as
            # given by 'grouped'.
            return op(0, expr)(seq(
                ifelse(lhs(t) = "none",
                       seq(thisproc(op(i, expr), marker), i in rhs(t)),
                       marker(op(0, expr)(seq(op(i, expr), i in rhs(t))))),
                    t in {entries(grouped, ':-pairs')}));
            
        elif type(expr, ':-anything' ^ ':-posint') then
            local finest_variable_sets := get_finest_variable_sets(op(1, expr));
            if `union`(op(variable_sets)) in variable_sets then
                return marker(expr);
            else
                return thisproc(op(1, expr), marker) ^ op(2, expr);
            end if;

        else
            return marker(expr);
        end if;    
    end proc;

local
    classifier := proc(i,
                       variable_sets_per_op,
                       maximal_nonsplittable,
                       $)
        for local mns in maximal_nonsplittable do
            if variable_sets_per_op[i] subset mns then
                return mns;
            end if;
        end do;

        return "none";
    end proc;

local
    get_finest_variable_sets := proc(expr, $)
    option cache;
        if type(expr, {':-`+`', ':-`*`'}) then
            local t;
            return `union`(seq(thisproc(t), t in expr));
        elif type(expr, ':-anything'^':-posint') then
            return thisproc(op(1, expr));
        else
            return {indets(expr, ':-name') minus {:-constants}};
        end if;
    end proc;
        
end module;

local
algebraic_expression_gen_diff :: static := proc(_self :: PowerSeriesObject,
                                                d :: nonnegint,
                                                $)
local ancestors := _self:-ancestors;
local n := numelems(ancestors:-variables);
local result := Array(1..0);

    # We should really never be called with d < max_derivative_computed + 1, because caching makes
    # sure the same thing is not recomputed. But if we *are* called in those circumstances, we need
    # to run the last iteration to compute the actual result. 
    for local dd from min(d, ancestors:-max_derivatives_computed + 1) to d do
        for local c in Iterator:-BoundedComposition([dd $ n], dd) do
            local i := ListTools:-SelectFirst(`>`, c, 0, 'output' = 'indices');
            local x := ancestors:-variables[i];
            # We compute derivatives[seq(c)] from derivatives[seq(cc)], where cc is obtained
            # from c by decreasing c[i] by 1, by differentiating w.r.t. x.
            c[i]--; # Now the value of c is cc as defined above...
            local parent := ancestors:-derivatives[seq(c)];
            # We want to compute:
            #
            #   diff(f(x), x $ (cc[i] + 1)) / (cc[i]+1)!,
            #
            # where f(x) is diff(_self:-algexpr, y$j, z$k) * y^j/j! * z^k/k!. We have:
            #
            #   parent = diff(f(x), x $ cc[i]) / cc[i]!.
            #
            # Then we need:
            #
            #   diff(parent, x) / (cc[i] + 1).
            local derivative := diff(parent, x) / (c[i] + 1);
            c[i]++; # ... and c is back to being itself.
            ancestors:-derivatives[seq(c)] := derivative;

            if dd = d then
                local j;
                result ,= eval(derivative, ancestors:-variables =~ 0)
                * mul(ancestors:-variables[j]^c[j], j=1..n);
            end if;
        end do;
        
        ancestors:-max_derivatives_computed++;
    end do;

    return add(result);
end proc;

local
algebraic_expression_gen_series :: static := proc(_self :: PowerSeriesObject,
                                                  d :: nonnegint,
                                                  $)
local ancestors := _self:-ancestors;
    if d < ancestors:-degrees_stored then
        return _self:-hpoly[d];
    end if;
local min_found_order := min(entries(ancestors:-found_order, ':-nolist'));
    if min_found_order = infinity then
        # The entry is infinity if there's no O(...) term, which means it's a polynomial. (Right?)
        # If that's the case for all variables, we have the full precision and any higher terms are
        # zero.
        return 0;
    end if;

local n := numelems(_self:-vars);
local target := ceil(max(min_found_order + 6, min_found_order * (n+1)/n, d));
local first_round := true;
    while min_found_order <= d do
        if first_round then
            first_round := false;
        else
            target += max(6, d - (min_found_order-1));
        end if;
    
        for local v in _self:-vars do
            ancestors:-used_order[v] := max(ancestors:-used_order[v],
                                            target + (ancestors:-used_order[v] - ancestors:-found_order[v]));
        end do;

        do_update_series(_self);

        min_found_order := min(entries(ancestors:-found_order, ':-nolist'));
    end do;

    return _self:-hpoly[d];
end proc;

local
do_update_series :: static := proc(_self :: PowerSeriesObject, $)
local ancestors := _self:-ancestors;
local expr := _self:-algexpr;
local minimum_order := infinity;

    userinfo(5, MultivariatePowerSeries, "calling do_update_series with used_order = %1",
             {indices(ancestors:-used_order, ':-pairs')});
    
    for local v in _self:-vars do
        expr := series(expr, v, ancestors:-used_order[v]);
        ancestors:-found_order[v] := get_order_from_series(expr, v);
        minimum_order := min(minimum_order, ancestors:-found_order[v]);
        expr := convert(expr, ':-polynom');
    end do;

    expr := expand(expr);
local split_per_degree := ListTools:-Classify(degree, convert(expr, ':-list', ':-`+`'), _self:-vars);
local loop_upper_bound := (if minimum_order = infinity then
                               max(indices(split_per_degree, ':-nolist'), upperbound(_self:-hpoly));
                           else
                               minimum_order - 1;
                           end if);

    if loop_upper_bound > ancestors:-degrees_stored then
        _self:-hpoly(loop_upper_bound + 1) := 0;
        for local d from ancestors:-degrees_stored + 1 to loop_upper_bound do
            if assigned(split_per_degree[d]) then
                _self:-hpoly[d] := add(split_per_degree[d]);
            else
                _self:-hpoly[d] := 0;
            end if;
            ancestors:-degrees_stored := d;
        end do;
        _self:-deg := loop_upper_bound;
    end if;
end proc;

local
get_order_from_series :: static := proc(expr,
                                        v :: name,
                                        $)
    if type(expr, ':-taylor') then
        if op(-2, expr) = ':-O'(1) then
            return op(-1, expr);
        else
            return infinity;
        end if;
        
    elif not type(expr, ':-series') then
        local terms := convert(expr, ':-list', ':-`+`');
        local bigO := select(type, terms, ':-O'(
            {':-constant', ':-identical'(v), ':-identical'(v)^':-nonnegint'}));

        if bigO = [] then
            # polynomial
            return infinity;

        elif numelems(bigO) = 1 then
            local bigO_argument := op(bigO[1]);

            if type(bigO_argument, ':-constant') then
                return 0;
            end if;

            bigO_argument := convert(bigO_argument, ':-list', ':-`^`');
            return bigO_argument[2];
        end if;
    end if;
    
    error "series did not produce a taylor expansion in %1", v;
end proc;
