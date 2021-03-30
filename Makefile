SRCS = PowerSeries/tst/addition.tst PowerSeries/tst/constant.tst PowerSeries/tst/deepCopy.tst PowerSeries/tst/fromPolynomial.tst PowerSeries/tst/fromProcedure.tst PowerSeries/tst/geometricSeries.tst PowerSeries/tst/identity.tst PowerSeries/tst/inverse.tst PowerSeries/tst/naryAdd.tst PowerSeries/tst/naryMul.tst PowerSeries/tst/one.tst PowerSeries/tst/product.tst PowerSeries/tst/sumOfAllMonomials.tst PowerSeries/tst/zero.tst
POWERSERIES_SRC = .
LIB_NAME = MultivariatePowerSeries

all: mla

mla: $(SRCS)
	rm -f $(LIB_NAME).mla  ../$(LIB_NAME).mla
	maple make_mla.mpl
	cp $(LIB_NAME).mla ../$(LIB_NAME).mla
	@(echo "Done!";)

test:
	$(foreach f, $(SRCS), maple $(f);)
	@(echo "Done!";)

clean:
	rm -f $(LIB_NAME).mla  ../$(LIB_NAME).mla

