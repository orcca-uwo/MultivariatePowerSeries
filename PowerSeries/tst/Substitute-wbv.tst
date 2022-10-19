#test 100

with(TestTools):
kernelopts(opaquemodules=false):
wbv := MultivariatePowerSeries:-PowerSeriesObject:-weighted_bounded_vectors:
kernelopts(opaquemodules=true):

# This doesn't really verify just as set; it also tests uniqueness *and* verifies that the sum and
# weight sum are correct at every step. For any *invalid* wbvs, it raises an error rather than
# returning false.
VerifyTools:-AddVerification(
    wbv_as_set = proc(w :: object,
                      s :: set,
                      weights :: list,
                      $)
    local ws := MutableSet();
    local v := w:-Values();
        
        do
            local i;
            if w:-Sum() <> add(v) then
                error "incorrect sum %1 for %2", w:-Sum(), v;
            elif w:-WeightSum() <> add(v[i] * weights[i], i = 1 .. numelems(v)) then
                error "incorrect sum %1 for %2 with weights %3", w:-WeightSum(), v, weights;
            end if;

            local vl := convert(v, 'list');
            if vl in ws then
                error "repeated entry: %1", v;
            end if;

            if not (vl in s) then
                return false;
            end if;
            ws:-insert(ws, vl);
        until not w:-Next();

        return andmap(`in`, s, ws);
    end proc);

pedestrian_wbv := proc(weights :: list(posint), wtsum :: nonnegint, $)
local i, j;
local n := numelems(weights);
local xs := [seq(cat(':-x', i), i=1..n)];
    return {value(foldl(%seq, [seq(xs[i], i = 1 .. n)],
                        seq(xs[i] = 0 .. (wtsum - add(weights[j] * xs[j], j = i+1 .. n))/weights[i],
                            i = 1 .. n)))};
end proc;

VerifyTools:-AddVerification(
    with_pedestrian_wbv = proc(o :: object,
                               params :: [list(posint), nonnegint],
                               $)
    local weights := params[1];
    local wtsum := params[2];
        return verify(o, pedestrian_wbv(weights, wtsum), 'wbv_as_set'(weights));
    end proc);
                             

# Length 0.
Try[verify,wbv_as_set([])]("0.0", Object(wbv, [], 0), {[]});
Try[verify,wbv_as_set([])]("0.1", Object(wbv, [], 5), {[]});

# Length 1.
Try[verify,wbv_as_set([1])]("1.0", Object(wbv, [1], 0), {[0]});
Try[verify,wbv_as_set([1])]("1.1", Object(wbv, [1], 5), {seq([i], i=0..5)});
Try[verify,wbv_as_set([3])]("1.2", Object(wbv, [3], 9), {seq([i], i=0..3)});
Try[verify,wbv_as_set([3])]("1.3", Object(wbv, [3], 10), {seq([i], i=0..3)});

# Length 2.
Try[verify,wbv_as_set([1,1])]("2.0", Object(wbv, [1,1], 0), {[0,0]});
Try[verify,wbv_as_set([1,1])]("2.1", Object(wbv, [1,1], 3), {
    [0, 0], [0, 1], [0, 2], [0, 3], [1, 0], [1, 1], [1, 2], [2, 0], [2, 1], [3, 0]});
Try[verify,wbv_as_set([1,3])]("2.2", Object(wbv, [1,3], 6), {
    [0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [2, 0],
    [2, 1], [3, 0], [3, 1], [4, 0], [5, 0], [6, 0]});
Try[verify,wbv_as_set([4,12])]("2.3", Object(wbv, [4,12], 24), {
    [0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [2, 0],
    [2, 1], [3, 0], [3, 1], [4, 0], [5, 0], [6, 0]});
Try[verify,wbv_as_set([1,3])]("2.4", Object(wbv, [1,3], 7), {
    [0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0],
    [2, 1], [3, 0], [3, 1], [4, 0], [4, 1], [5, 0], [6, 0], [7, 0]});

# Length 3.
Try[verify,wbv_as_set([1,1,1])]("3.0", Object(wbv, [1,1,1], 0), {[0,0,0]});
Try[verify,wbv_as_set([1,3,5])]("3.1", Object(wbv, [1,3,5], 0), {[0,0,0]});
Try[verify,wbv_as_set([1,3,5])]("3.2", Object(wbv, [1,3,5], 6), {
    [0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 2, 0], [1, 0, 0], [1, 0, 1], [1, 1, 0],
    [2, 0, 0], [2, 1, 0], [3, 0, 0], [3, 1, 0], [4, 0, 0], [5, 0, 0], [6, 0, 0]});
Try[verify,wbv_as_set([3,3,3])]("3.3", Object(wbv, [3,3,3], 4), {
    [0, 0, 0], [0, 0, 1], [0, 1, 0], [1, 0, 0]});

# Now test all of these things vs pedestrian_wbv, so that we can trust it at least somewhat for the
# final group of tests:
# Length 0.
Try("4.0.0", pedestrian_wbv([], 0), {[]});
Try("4.0.1", pedestrian_wbv([], 5), {[]});
# Length 1.
Try("4.1.0", pedestrian_wbv([1], 0), {[0]});
Try("4.1.1", pedestrian_wbv([1], 5), {seq([i], i=0..5)});
Try("4.1.2", pedestrian_wbv([3], 9), {seq([i], i=0..3)});
Try("4.1.3", pedestrian_wbv([3], 10), {seq([i], i=0..3)});
# Length 2.
Try("4.2.0", pedestrian_wbv([1,1], 0), {[0,0]});
Try("4.2.1", pedestrian_wbv([1,1], 3), {
    [0, 0], [0, 1], [0, 2], [0, 3], [1, 0], [1, 1], [1, 2], [2, 0], [2, 1], [3, 0]});
Try("4.2.2", pedestrian_wbv([1,3], 6), {
    [0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [2, 0],
    [2, 1], [3, 0], [3, 1], [4, 0], [5, 0], [6, 0]});
Try("4.2.3", pedestrian_wbv([4,12], 24), {
    [0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [2, 0],
    [2, 1], [3, 0], [3, 1], [4, 0], [5, 0], [6, 0]});
Try("4.2.4", pedestrian_wbv([1,3], 7), {
    [0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0],
    [2, 1], [3, 0], [3, 1], [4, 0], [4, 1], [5, 0], [6, 0], [7, 0]});
# Length 3.
Try("4.3.0", pedestrian_wbv([1,1,1], 0), {[0,0,0]});
Try("4.3.1", pedestrian_wbv([1,3,5], 0), {[0,0,0]});
Try("4.3.2", pedestrian_wbv([1,3,5], 6), {
    [0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 2, 0], [1, 0, 0], [1, 0, 1], [1, 1, 0],
    [2, 0, 0], [2, 1, 0], [3, 0, 0], [3, 1, 0], [4, 0, 0], [5, 0, 0], [6, 0, 0]});
Try("4.3.3", pedestrian_wbv([3,3,3], 4), {
    [0, 0, 0], [0, 0, 1], [0, 1, 0], [1, 0, 0]});

# Longer lengths.
Try[verify,with_pedestrian_wbv]("5.0", Object(wbv, [1,2,3,4,5], 10), [[1,2,3,4,5], 10]);
Try[verify,with_pedestrian_wbv]("5.1", Object(wbv, [1,1,2,3,5,8], 13), [[1,1,2,3,5,8], 13]);
Try[verify,with_pedestrian_wbv]("5.2", Object(wbv, [1$10], 10), [[1$10], 10]);

#end test
