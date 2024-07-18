######################
# ## Internal Routines 
######################

# ensure up to degree 'deg' is generated in hpoly,
# if not, then this routine will generate the missing parts
# and update self:-hpoly and self:-deg using self:-gen. 
# Note deg = -1 is the debugging signal.  
local
    ensure_degree ::static := proc (self :: PowerSeriesObject,
                                    deg :: {nonnegint, -1},
                                    $)
        if self:-deg < deg then
            if upperbound(self:-hpoly) < deg then
                self:-hpoly(deg+1) := 0;
            end if;
            for local i from self:-deg + 1 to deg do
                self:-deg := i; 
                self:-hpoly[i] := self:-gen(self, i);
            end do;
        end if;
        return self;
    end proc;

# convert poly to hpoly::Array  
local 
    convert_hpoly_from_poly ::static := proc(p :: polynom,
                                             {output :: identical(table, Array) := Array},
                                             $)
        local terms := p; # p is already expanded!
        local ld := degree(terms);
        terms := convert(terms, ':-list', ':-`+`');
        if not type(terms, ':-list'({polynom})) then
            error "unsupported term: %1", remove(type, terms, {polynom})[1];
        elif terms = [0] then
            terms := table();
        else
            terms := ListTools:-Classify(degree, terms);
        end if;
        
        if output = ':-table' then
            return map(add, eval(terms, 1));
        else
            return Array(0 .. ld, map(eq -> lhs(eq) = add(rhs(eq)),  {indices(terms, 'pairs')}));
        end if;
    end proc;
    
# get d-th homogeneous part of polynomial p
local 
    get_homog_part_of_poly ::static := proc(p :: polynom,
                                            d :: nonnegint,
                                            $)
        local dummy;
        local terms := AUTO_EXPAND(dummy, p);
        local ld := degree(terms);
        if ld < d then return 0; end;
        terms := convert(terms, ':-list', ':-`+`');
        if not type(terms, ':-list'({polynom})) then
            error "unsupported term: %1", remove(type, terms, {polynom})[1];
        end if;
        return add(select(t -> degree(t) = d, terms));
    end proc;

###############
# ## Generators 
###############

# poly_gen
local 
    poly_gen ::static := proc(self :: PowerSeriesObject, 
                        d :: nonnegint,
                        $)
        return ifelse(d <= self:-deg, self:-hpoly[d], 0)
    end proc;

# proc_gen
local  
    proc_gen ::static := proc(self :: PowerSeriesObject, 
                d :: nonnegint,
                $)
        local p_d := self:-ancestors:-P(d);
        if self:-ancestors:-E then 
            local dummy;
            p_d := AUTO_EXPAND(dummy, p_d);
        end if;
        if self:-ancestors:-V and p_d <> 0 then
            local monomials := convert(p_d, list, ':-`+`');
            local vars := self:-vars;
            monomials := remove(t -> type(t, ':-polynom'(':-constant', vars))
                                and degree(t, vars) = d, monomials);
            if monomials <> [] then
                error "the generator procedure was expected to return a polynomial in %1 of total degree "
                "%2, but returned %3: %4",
                vars, d, ifelse(numelems(monomials) = 1, "this term", "these terms"), add(monomials);
            end if;
        end if;
        return p_d;
    end proc;

# add_gen
local 
    add_gen ::static := proc(self :: PowerSeriesObject,  
                            d :: nonnegint,
                            $) 
        return HOMOGENEOUS_PART(self:-ancestors:-A, d) + HOMOGENEOUS_PART(self:-ancestors:-B, d);
    end proc;

# nary_add_gen
local 
    nary_add_gen ::static := proc(self :: PowerSeriesObject, 
                                d ::nonnegint, 
                                $)
        local t;
        return add(HOMOGENEOUS_PART(t, d), t in self:-ancestors:-T);
    end proc;

# nary_add_with_coeffs_gen
local 
    nary_add_with_coeffs_gen ::static := proc(self :: PowerSeriesObject, 
                                            d ::nonnegint, 
                                            $)
        local p, dummy;
        return add(AUTO_EXPAND(dummy, HOMOGENEOUS_PART(p[1], d) * p[2]), p in self:-ancestors:-P);
    end proc;

# nary_add_with_poly_coeffs_gen
local 
    nary_add_with_poly_coeffs_gen ::static := proc(self :: PowerSeriesObject, 
                                                d ::nonnegint, 
                                                $)
    local pp;
        return add(HOMOGENEOUS_PART(pp, d), pp in self:-ancestors:-S);
    end proc;

# sub_gen
local 
    sub_gen ::static := proc(self :: PowerSeriesObject, 
                            d :: nonnegint,
                            $)
            return HOMOGENEOUS_PART(self:-ancestors:-A, d) - HOMOGENEOUS_PART(self:-ancestors:-B, d);
    end proc;

# neg_gen
local 
    neg_gen ::static := proc(self :: PowerSeriesObject, 
                    d :: nonnegint,
                    $)
        return - HOMOGENEOUS_PART(self:-ancestors:-A, d);
    end proc;

# mul_gen
local 
    mul_gen ::static := proc(self :: PowerSeriesObject, 
                    d :: nonnegint,
                    $)
        local dummy, i;
        return AUTO_EXPAND(dummy, add( ifelse( HOMOGENEOUS_PART(self:-ancestors:-A, i) = 0, 0, HOMOGENEOUS_PART(self:-ancestors:-A, i) * HOMOGENEOUS_PART(self:-ancestors:-B, d-i)), i = 0 .. d));
    end proc;

# Note: using Iterator makes mul computation of unbalanced inputs very slow because of intesively usage of Maple expand    
# local 
#     nary_mul_gen ::static := proc(self :: PowerSeriesObject, 
#                                 d :: nonnegint, 
#                                 $)
#         local n := numelems(self:-ancestors:-F);
#         # local i;
#         # local result := 0;
#         local dummy;
#         local result := Array(1..0);
#         local homoparts := Array(1 .. n);
#         for local degree_array in Iterator:-BoundedComposition([d $ n], d) do
#             # ArrayTools:-Append(result, AUTO_EXPAND(dummy, mul(HOMOGENEOUS_PART(self:-ancestors:-F[i], degree_array[i]), i = 1 .. n))); # speed-up ~1.4 instead of doing partial adds..
#             for local j from 1 to n do 
#                 homoparts[j] := HOMOGENEOUS_PART(self:-ancestors:-F[j], degree_array[j]);
#                 next degree_array if homoparts[j] = 0;
#             end do;
#             ArrayTools:-Append(result, AUTO_EXPAND(dummy, mul(homoparts)));
#         end do;
#         return add(result);
#     end proc;

# beq_gen
local
    beq_gen ::static := proc(self :: PowerSeriesObject,
                    d :: nonnegint,
                    $)
        local s := HOMOGENEOUS_PART(self:-ancestors:-A, d);
        self:-ensure_degree(self, d-1);
        local i, dummy;
        s -= add(AUTO_EXPAND(dummy, HOMOGENEOUS_PART(self:-ancestors:-B, i) * self:-hpoly[d-i]), i = 1 .. d);
        return AUTO_EXPAND(dummy, normal(s / self:-ancestors:-B:-hpoly[0]));
    end proc;

# To convert a PuSO to a PSO    
# ConvertToPSO_gen
local ConvertToPSO_gen ::static := proc(_self :: PowerSeriesObject, 
                                        d :: nonnegint, $) 
    local x, result, terms, terms_as_products, terms_as_power_products,
          total_degrees, new_entry, lower_bnd, degrees_as_list;
    local ancestors := _self:-ancestors;
    local A := ancestors:-A;

    if ancestors:-is_error then
        error "the Puiseux series used to create %1 is not"
              " a power series", _self;
    end if;

    local bound := A:-AbsoluteDegreeBound(A, d);
    local pso := A:-GetPowerSeries(A);
    local anc_size := upperbound(ancestors:-cache_table);

    if bound>anc_size then
        ancestors:-cache_table(bound+1) := 0;  
    end if;

    lower_bnd := anc_size+1;

    for local i from lower_bnd to bound do
        result := pso:-HomogeneousPart(pso, i);
        result := subs(A:-ChangeOfVariables(A), result);
        
        # Create a list of the monomials.
        terms := convert(result, ':-list', ':-`+`');
        
        # "Multiply" by the monomial given by e
        terms := map(':-`*`', terms, A:-GetMonomial(A));
        
        # Write each element in 'terms' as a list of the factors.
        terms_as_products := map(convert, terms, ':-list', ':-`*`');
        
        # Write each factor in the sub-lists of 'terms_as_products' as a 
        #list of the base and exponent.
        terms_as_power_products := map2(map, convert, terms_as_products, ':-list', ':-`^`');
        
        # Keep in 'terms_as_power_products' only the sub-sub-lists that correspond to variables 
        #(i.e. remove the coefficient terms).
        terms_as_power_products := map2( select, type, terms_as_power_products, [':-name',':-rational'] );

        degrees_as_list := map[3]( map[2], op, 2, terms_as_power_products );
        # Determine the total absolute degrees of each monomial in 'terms'.
        total_degrees := map( add, degrees_as_list);

        if not type(degrees_as_list, ':-list'(':-list'(':-nonnegint'))) then
            ancestors:-is_error := true;
            _self:-hpoly := Array(0..0, [undefined]);
            error "the Puiseux series used to create %1 is not"
                  " a power series", _self;
        end if;

        new_entry := table(':-sparse');
        for local j to numelems(terms) do
          new_entry[total_degrees[j]] += terms[j];
        end do;

        ancestors:-cache_table[i] := eval(new_entry,1);
    end do;

    local output := 0;
    for local i from 0 to bound do
        output += ancestors:-cache_table[i][d];
    end do;

    return output;
end proc;

# Generator for the function univariate_product_inverse_monomial.
local prod_inv_mon_gen ::static := proc(_self :: PowerSeriesObject, 
                                         d :: nonnegint, $) 
 
    local A := _self:-ancestors:-A;
    local q := _self:-ancestors:-q;
    local deg := degree(q);

    return HOMOGENEOUS_PART(A, d+deg)/q; 
end proc;

local force_precision ::static := proc(self :: PowerSeriesObject,
                    d :: nonnegint,
                    $)
    local real_precision := numelems(self:-hpoly)-1;
    if d > real_precision then 
        error "forced degree should be no larger than the real "
                "precision, %1, but received %2", real_precision, d;
    end;
    local old_precision := self:-deg;
    self:-deg := d;

    return old_precision;
end proc;

local real_precision ::static := proc(self :: PowerSeriesObject,
                    $)
    return numelems(self:-hpoly)-1;
end proc;