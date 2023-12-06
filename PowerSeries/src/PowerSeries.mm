module PowerSeriesObject()
option object;

    # PowerSeriesObject Module Attributes
    local 
        # hpoly : To store the array of homogeneous parts started at index = 0
        #         This is an Maple Array to acheive performance in  updating, reshaping, copying, and indexing. 
        hpoly :: Array,
        # deg : To store the max degree of computed homogeneous parts 
        #       It basically equals to upperbound(hpoly), unless we know the required precision 
        #       in advanced and do some optimization (e.g. ensure_degree) 
        deg :: nonnegint,
        # gen : Every formal power series must have a generator to compute the homogeneous part
        #       of 0 <= index < infinity. see generator local constructors in internal_routines.mm
        gen :: procedure,
        # vars : set containing all variable names occurring in this power series.
        vars :: set(name),
        # ancestors : data used by a generator 
        # Note we add this local as self:-gen cannot use local variables where it defines (surprisingly!),
        # using ancestors is an alternative to fix this issue, instead of using the inheritance 
        # mechanism in which we faced a bug (MPL-23960). 
        ancestors :: record,
        # algexpr : Analytical Expression is in-fact a helper attribute used mainly in Display and ModulePrint 
        #           and presents the analytic expression that the series tends towards (in a formal
        #           sense at least)
        algexpr :: algebraic,
        # dstyle : a display style for a specific power series  
        dstyle :: PSDS_TYPE; 



# return an in-place subArray ps:-hpoly[0 .. d]. Note d <= ps:-deg 
$define UP_TO_DEGREE_ARRAY(p, d) \
    ifelse(p:-deg = d, p:-hpoly, ArrayTools:-Alias(p:-hpoly, [0 .. d]))

# return homogeneous part of deg = d 
# Note if d > ps:-deg this macro will generate (d - ps:-deg) homogeneuos parts 
$define HOMOGENEOUS_PART(p, d) \
    ifelse(p:-deg >= d, p:-hpoly[d], p:-ensure_degree(p, d):-hpoly[d])

# return Expanded of x wehere d is a dummy variable for the sake of efficiency 
$define AUTO_EXPAND(d, x) \
    ifelse(type((d := x), ':-polynom'(':-complex'(':-numeric'))), expand(d), evala(':-Expand'(d)))
 
# ModuleApply: to make PowerSeriesObject an object factory
local 
    ModuleApply ::static := proc()
        Object(PowerSeriesObject, _passed);
    end proc;

# ModuleCopy: is called when a new object of the same class is created via a call to Object
local
    ModuleCopy ::static := proc(new :: PowerSeriesObject,
                                old :: PowerSeriesObject,
                                hpoly :: Array,
                                deg :: nonnegint,
                                gen :: procedure,
                                vars :: set(name),
                                ancestors :: {list({name, string} = anything), record} := [],
                                algexpr :: algebraic := undefined,
                                dstyle :: PSDS_TYPE := [],
                                {initializeconstantterm::truefalse :=false},
                                $)
        if _npassed > 2 then
            ASSERT(lowerbound(hpoly) = 0, "array with wrong lowerbound");
            ASSERT(upperbound(hpoly) = deg, "wrong degree specified");
            new:-hpoly := hpoly;
            new:-deg := deg;
            new:-gen := gen;
            new:-vars := vars;
            if type(ancestors, record) then
                new:-ancestors := ancestors;
            else
                new:-ancestors := Record[':-packed'](op(ancestors));
            end if;
            new:-dstyle := dstyle;
            new:-algexpr := algexpr;
        elif type(old:-hpoly, 'Array') then
            new:-hpoly := copy(old:-hpoly);
            new:-deg := old:-deg;
            new:-gen := old:-gen;
            new:-ancestors := old:-ancestors;
            new:-dstyle := old:-dstyle;
            new:-algexpr := old:-algexpr;
        else
            error "you cannot copy the original PowerSeriesObject object";
        end if;

        if initializeconstantterm then
            new:-hpoly[0] := new:-gen(new, 0);
        end if;
    end proc;
    
local 
    defaultDisplayStyle ::static := [];

# SetDefaultDisplayStyle
# To set the local global variable defaultDisplayStyle 
export 
    SetDefaultDisplayStyle ::static := proc(s :: PSDS_TYPE, $)
        defaultDisplayStyle := s;
    end proc;
# GetDefaultDisplayStyle
export 
    GetDefaultDisplayStyle ::static := proc($)
        return defaultDisplayStyle;
    end proc;

# SetDisplayStyle
# to set the dstyle of the input object
export    
    SetDisplayStyle ::static := proc(obj :: PowerSeriesObject, s :: PSDS_TYPE, $)
        obj:-dstyle := s; 
    end proc;

local
  typeset_sum_of_terms :: static := proc(l :: list,
                                         typeset_tail := undefined,
                                         $)
      uses T = Typesetting;
      local result; # expression sequence of terms to be typeset
      if l = [] then
          if typeset_tail = undefined then
              result := T:-mn("0");
          else
              result := typeset_tail;
          end if;
      else
          result := Array(1..1, [apply(T:-Typeset, l[1])]);
          for local term in l[2..] do
              # if type(term, ':-negative') or (type(term, ':-`*`') and membertype(':-negative', [op(term)])) then                
              if `tools/sign`(term) = -1 then
                  result ,= T:-mo("&minus;");
                  # just to be safe...
                  if type(-term, ':-`+`') then
                      result ,= T:-mfenced(apply(T:-Typeset, -term));
                  else
                      result ,= apply(T:-Typeset, -term);
                  end if;
              else
                  result ,= T:-mo("&plus;");
                  result ,= apply(T:-Typeset, term);
              end if;
          end do;
          if typeset_tail <> undefined then
              result ,= T:-mo("&plus;");
              result ,= typeset_tail;
          end if;
          result := seq(result);
      end if;
      return result;
  end proc;

local
  sprintf_sum_of_terms :: static := proc(l :: list, $)
      if l = [] then
          return "0";
      end if;

  local result := Array(1..1, [sprintf("%a", l[1])]);
      for local term in l[2..] do
          if `tools/sign`(term) = -1 then
              result ,= ':-`-`';
              # just to be safe...
              if type(-term, ':-`+`') then
                  result ,= sprintf("(%a)", -term);
              else
                  result ,= sprintf("%a", -term);
              end if;
          else
              result ,= ':-`+`';
              result ,= sprintf("%a", term);
          end if;
      end do;

      return StringTools:-Join(convert(result, ':-list'), " ");
  end proc;
    
# Display
# to disply the input object, 
# user_dstyle : will overwrite the default and self:-dstyle style 
# output : (doesn't need to be documented) used for internal purposes
export 
    Display ::static := proc(self :: PowerSeriesObject, 
                             user_dstyle :: PSDS_TYPE := [],
                             output :: identical("full", "typeset", "fullstring", "string") := "full",
                             {updateterms :: name := NULL},
                             {substitutions :: {set,list}(equation) := {}},
                             {objectname :: string := "PowerSeries"},
                             {cofactor :: algebraic := 1},
                             $)
        uses T = Typesetting;
        local output_string := member(output, {"string", "fullstring"});
        local mystyle;
        if user_dstyle <> [] then
            mystyle := user_dstyle;
        elif self:-dstyle <> [] then
            mystyle := self:-dstyle;
        else
            mystyle := defaultDisplayStyle;
        end if;

        local terms := eval(':-maxterms', [op(mystyle), ':-maxterms'=50]);
        local min_prec := min(self:-deg, eval(':-precision', [op(mystyle), ':-precision'=infinity]));
        local result := Array(0 .. min_prec, ':-fill' = NULL);
        local nterms := 0; 
        local tail; 
        local partial_degree := false;

        local d;
        for d from 0 to min_prec while nterms < terms do
            local h := subs(substitutions, self:-hpoly[d]);
            if h = 0 then next; end if;
            local hl := convert(h, list, ':-`+`');
            hl := map(:-`*`, hl, cofactor);
            if nterms + numelems(hl) > terms then
                tail := hl[1..terms-nterms];
                nterms := terms;
                partial_degree := true;
                break;
            else 
                result[d] := op(hl);
                nterms += numelems(hl);
            end if;
        end do;

        if updateterms <> NULL then
            updateterms -= nterms;
        end if;

        local result_is_known_complete := (not partial_degree) and self:-algexpr <> undefined
        and type(self:-algexpr, ':-polynom') and d > degree(self:-algexpr);
        result := convert(result, ':-list');

        if partial_degree then
            if add(result) = 0 then
                if output_string then
                    result := cat(sprintf_sum_of_terms(tail), " + ...");
                else
                    result := typeset_sum_of_terms(tail), T:-mo("&plus;"), T:-mn("&hellip;");
                end if;
            else
                if output_string then
                    result := cat(sprintf_sum_of_terms(result), " + (",
                                  sprintf_sum_of_terms(tail), " + ...) + ...");
                else
                    tail := T:-mfenced(T:-mrow(typeset_sum_of_terms(tail), T:-mo("&plus;"), T:-mn("&hellip;")));
                    result := typeset_sum_of_terms(result, tail), T:-mo("&plus;"), T:-mn("&hellip;");
                end if;
            end if;

        elif result_is_known_complete then
            if output_string then
                result := sprintf_sum_of_terms(result);
            else
                result := typeset_sum_of_terms(result);
            end if;
            
        else
            if output_string then
                result := cat(sprintf_sum_of_terms(result), " + ...");
            else
                result := typeset_sum_of_terms(result), T:-mo("&plus;"), T:-mn("&hellip;");
            end if;
        end if;

        if not member(output, {"typeset", "string"}) then
            local alg := subs(substitutions, self:-algexpr);
            if alg = undefined or (type(alg, polynom) and result_is_known_complete) then
                if output_string then
                    result := sprintf("< %s: %s >", objectname, result);
                else
                    result := T:-mfenced(T:-mrow( T:-mn(sprintf("%s:   ", objectname)), result), ':-open' = "&lsqb;", ':-close' = "&rsqb;");
                end if;
            else
                alg *= cofactor;
                if output_string then
                    local algstr := sprintf("%a", alg);
                    if length(algstr) > 100 then
                        algstr := cat(algstr[.. 10], " ... ", algstr[-10 ..]);
                    end if;
                    result := sprintf("< %s of %s : %s >", objectname, algstr, result);
                else
                    local algt := apply(T:-Typeset, alg);
                    if length(algt) > 500 then
                        algt := elide_typeset(algt, 500);
                    end if;
                    result := T:-mfenced(T:-mrow( T:-mn(sprintf("%s of ", objectname)), algt, T:-mn(" : "), result), ':-open' = "&lsqb;", ':-close' = "&rsqb;");
                end if;
            end if;
        end if;

        return result;
    end proc;

local
    elide_typeset :: static := proc(t,
                                    max_length :: integer,
                                    $)
    uses T = Typesetting;
    local length_left; # assigned in a branch in the if statement

        if length(t) < max_length or nops(t) = 0 then
            return t;
            
        elif (type(t, ':-anyfunc'(':-string')) and length(t) > max_length)
        or not type(t, ':-specfunc'({
            T:-mn, T:-mi, T:-mo, T:-mio, T:-mtext, T:-mspace, T:-ms, T:-mglyph, T:-mrow, T:-mfrac,
            T:-msqrt, T:-mroot, T:-mfenced, T:-msub, T:-msup, T:-msubsup, T:-munder, T:-mover,
            T:-munderover}))
        or (length_left := max_length - length(op(0, t)) - 10) <= 5 then
            return T:-mn("&hellip;");

        elif nops(t) = 1 then
            return applyop(thisproc, 1, t, length_left);
        end if;

        # Now 5 < length_left < max_length <= length(t), nops(t) > 1, and t is a specfunc of one of
        # the listed types.
    local n := nops(t);
    local lengths := map(length, [op(t)]);

        # Is there one entry that takes it over the threshold?
        if add(lengths) - max(lengths) < length_left - 5 then
            local max_index := n+1 - max[':-index'](ListTools:-Reverse(lengths));
            return applyop(thisproc, max_index, t, length_left - 5 - (add(lengths) - max(lengths)));
        end if;

        # We keep adding operands, one by one, to what we will return; starting with the first
        # operand, then the last one, then the ones in between. If we include an operand and an
        # operator follows it (precedes it if we're looking at the last entry), we include it, too.
    local include_start := 0;
    local include_end := 0;
        do
            local this_op, next_op;
            if include_start = 0 then
                this_op, next_op := 1, 2;
            elif include_end = 0 then
                this_op, next_op := n, n-1;
            else
                this_op, next_op := include_start + 1, include_start + 2;
            end if;

            if lengths[this_op] > length_left - 5 then
                # Recurse on this_op, elide what's left.
                if type(op(next_op, t), ':-specfunc'(T:-mo)) then
                    # Include next_op, too.
                    if this_op = n then
                        return op(0, t)(
                            op(1 .. include_start, t),
                            ifelse(include_start + 2 < n, T:-mn("&hellip;"), NULL),
                            op(next_op, t),
                            thisproc(op(this_op, t), length_left - 5 - lengths[next_op]));
                    else
                        return op(0, t)(
                            op(1 .. include_start, t),
                            thisproc(op(this_op, t), length_left - 5 - lengths[next_op]),
                            op(next_op, t),
                            ifelse(include_start + include_end + 2 < n, T:-mn("&hellip;"), NULL),
                            op(n + 1 - include_end .. n, t));
                    end if;
                else
                    # Ignore next_op.
                    if this_op = n then
                        return op(0, t)(
                            op(1 .. include_start, t),
                            ifelse(include_start + 1 < n, T:-mn("&hellip;"), NULL),
                            thisproc(op(this_op, t), length_left - 5));
                    else
                        return op(0, t)(
                            op(1 .. include_start, t),
                            thisproc(op(this_op, t), length_left - 5),
                            ifelse(include_start + include_end + 1 < n, T:-mn("&hellip;"), NULL),
                            op(n + 1 - include_end .. n, t));
                    end if;
                end if;
            end if;

            if type(op(next_op, t), ':-specfunc'(T:-mo)) then
                # Include this_op and next_op
                length_left -= lengths[this_op] + lengths[next_op];
                if this_op = n then
                    include_end := 2;
                else
                    include_start += 2;
                end if;
            else
                # Include this_op; ignore next_op
                length_left -= lengths[this_op];
                if this_op = n then
                    include_end := 1;
                else
                    include_start += 1;
                end if;
            end if;
        end do;
    end proc;

# ModulePrint: print self :: PowerSeriesObject 
local
    ModulePrint :: static := proc(self :: PowerSeriesObject, $)
        if IsWorksheetInterface() then
            return Display(self);
        else
            return convert(Display(self, [], "fullstring"), ':-name');
        end if;
    end proc;

# Module Deconstruct
local
  ModuleDeconstruct :: static := proc(self :: PowerSeriesObject, $)
    return 'Object'('PowerSeriesObject', self:-hpoly, self:-deg, self:-gen, self:-vars, self:-ancestors);
  end proc;

$include "MultivariatePowerSeries/PowerSeries/src/internal_routines.mm"
$include "MultivariatePowerSeries/PowerSeries/src/basic_routines.mm"
$include "MultivariatePowerSeries/PowerSeries/src/basic_arithmetic.mm"
$include "MultivariatePowerSeries/PowerSeries/src/FactorOutMonomial.mm"
$include "MultivariatePowerSeries/PowerSeries/src/FromAlgebraicExpression.mm"
$include "MultivariatePowerSeries/PowerSeries/src/Substitute.mm"
$include "MultivariatePowerSeries/PowerSeries/src/TaylorShift.mm"

$undef UP_TO_DEGREE_ARRAY
$undef HOMOGENEOUS_PART
$undef AUTO_EXPAND

end module:
