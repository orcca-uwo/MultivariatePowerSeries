# Exponentiation
# to compute _self^n 
# n : nonnegint 
export Exponentiate ::static := proc(_self :: PuiseuxSeriesObject, 
                                n :: nonnegint, 
                                $)
        if n = 0 then 
            return Object(PuiseuxSeriesObject, 1); 
        elif n = 1 then 
            return _self;
        end if;        
        local p := _self;
        local m := abs(n);
        local res := Array(0 .. -1);
        while m > 0 do 
            if irem(m,2) = 1 then 
                res ,= p;
            end if;
            p := p:-BinaryMultiply(p, p);
            m := iquo(m,2);
        end do;
        return p:-NaryMultiply(seq(res));
    end proc;