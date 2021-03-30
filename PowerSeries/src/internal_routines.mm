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
            self:-hpoly(deg+1) := 0; 
            for local i from self:-deg + 1 to deg do
                self:-hpoly[i] := self:-gen(self, i);
                self:-deg := i; 
            end do;
        end if;
        return self;
    end proc;

# convert poly to hpoly::Array  
local 
    convert_hpoly_from_poly ::static := proc(p :: polynom, $)
        local terms := p; # p is already expanded!
        local ld := degree(terms);
        terms := convert(terms, ':-list', ':-`+`');
        if not type(terms, ':-list'({polynom})) then
            error "unsupported term: %1", remove(type, terms, {polynom})[1];
        end if;
        terms := ListTools:-Classify(degree, terms);
        return Array(0 .. ld, map(eq -> lhs(eq) = add(rhs(eq)),  {indices(terms, 'pairs')}));
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
        return HomogeneousPart(self:-ancestors:-A, d) + HomogeneousPart(self:-ancestors:-B, d);
    end proc;

# nary_add_gen
local 
    nary_add_gen ::static := proc(self :: PowerSeriesObject, 
                                d ::nonnegint, 
                                $)
        local t;
        return add(HomogeneousPart(t, d), t in self:-ancestors:-T);
    end proc;

# nary_add_with_coeffs_gen
local 
    nary_add_with_coeffs_gen ::static := proc(self :: PowerSeriesObject, 
                                            d ::nonnegint, 
                                            $)
        local p, dummy;
        return add(AUTO_EXPAND(dummy, HomogeneousPart(p[1], d) * p[2]), p in self:-ancestors:-P);
    end proc;

# nary_add_with_poly_coeffs_gen
local 
    nary_add_with_poly_coeffs_gen ::static := proc(self :: PowerSeriesObject, 
                                                d ::nonnegint, 
                                                $)
    local pp;
        return add(HomogeneousPart(pp, d), pp in self:-ancestors:-S);
    end proc;

# sub_gen
local 
    sub_gen ::static := proc(self :: PowerSeriesObject, 
                            d :: nonnegint,
                            $)
            return HomogeneousPart(self:-ancestors:-A, d) - HomogeneousPart(self:-ancestors:-B, d);
    end proc;

# neg_gen
local 
    neg_gen ::static := proc(self :: PowerSeriesObject, 
                    d :: nonnegint,
                    $)
        return - HomogeneousPart(self:-ancestors:-A, d);
    end proc;

# mul_gen
local 
    mul_gen ::static := proc(self :: PowerSeriesObject, 
                    d :: nonnegint,
                    $)
        local dummy, i;
        return AUTO_EXPAND(dummy, add( ifelse( HomogeneousPart(self:-ancestors:-A, i) = 0, 0, HomogeneousPart(self:-ancestors:-A, i) * HomogeneousPart(self:-ancestors:-B, d-i)), i = 0 .. d));
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
#             # ArrayTools:-Append(result, AUTO_EXPAND(dummy, mul(HomogeneousPart(self:-ancestors:-F[i], degree_array[i]), i = 1 .. n))); # speed-up ~1.4 instead of doing partial adds..
#             for local j from 1 to n do 
#                 homoparts[j] := HomogeneousPart(self:-ancestors:-F[j], degree_array[j]);
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
        local s := HomogeneousPart(self:-ancestors:-A, d);
        self:-ensure_degree(self, d-1);
        local i, dummy;
        s -= add(AUTO_EXPAND(dummy, HomogeneousPart(self:-ancestors:-B, i) * self:-hpoly[d-i]), i = 1 .. d);
        return AUTO_EXPAND(dummy, normal(s / self:-ancestors:-B:-hpoly[0]));
    end proc;
