SRCS = MultivariatePowerSeries/src/MultivariatePowerSeries.mpl
MLAFILE = MultivariatePowerSeries.mla
ASSERTLEVEL = 2

ALLTESTFILES=PuiseuxSeries/tst/addition.tst PuiseuxSeries/tst/convert_to_pso.tst PuiseuxSeries/tst/exports.tst PuiseuxSeries/tst/inverse.tst PuiseuxSeries/tst/longerExample1.tst PuiseuxSeries/tst/longerExample2.tst PuiseuxSeries/tst/makeCompatible.tst PuiseuxSeries/tst/multiplication.tst PuiseuxSeries/tst/objectDefinition.tst PuiseuxSeries/tst/randomtest_MakeRaysCompatible.tst PuiseuxSeries/tst/randomtest_objectDefinition.tst PuiseuxSeries/tst/trigonometricIdentity1.tst PuiseuxSeries/tst/trigonometricIdentity2.tst PuiseuxSeries/tst/truncate.tst PowerSeries/tst/addition.tst PowerSeries/tst/constant.tst PowerSeries/tst/deepCopy.tst PowerSeries/tst/FactorOutMonomial.tst PowerSeries/tst/fromPolynomial.tst PowerSeries/tst/fromProcedure.tst PowerSeries/tst/geometricSeries.tst PowerSeries/tst/identity.tst PowerSeries/tst/inverse.tst PowerSeries/tst/naryAdd.tst PowerSeries/tst/naryMul.tst PowerSeries/tst/one.tst PowerSeries/tst/product.tst PowerSeries/tst/Substitute.tst PowerSeries/tst/Substitute-wbv.tst PowerSeries/tst/sumOfAllMonomials.tst PowerSeries/tst/TaylorShift.tst PowerSeries/tst/zero.tst UPoPS/tst/addition.tst UPoPS/tst/constant.tst UPoPS/tst/fromPolynomial.tst UPoPS/tst/hensel_factorization.tst UPoPS/tst/one.tst UPoPS/tst/product.tst UPoPS/tst/Puiseux_theorem.tst UPoPS/tst/taylorShift.tst UPoPS/tst/UPoPSOverPuSO.tst UPoPS/tst/wp.tst UPoPS/tst/zero.tst

SOMETESTFILES = PowerSeries/tst/addition.tst PowerSeries/tst/constant.tst PowerSeries/tst/deepCopy.tst PowerSeries/tst/FactorOutMonomial.tst PowerSeries/tst/fromPolynomial.tst PowerSeries/tst/fromProcedure.tst PowerSeries/tst/geometricSeries.tst PowerSeries/tst/identity.tst PowerSeries/tst/inverse.tst PuiseuxSeries/tst/addition.tst

ALLTESTARGS = $(ALLTESTFILES:%=%.test)

SANITYTARGS = $(SOMETESTFILES:%=%.test)

%.test:
	maple -A 2 -B -b MultivariatePowerSeres.mla $*

.PHONY: clean mint test mla

all: mla

mint: 
	@$(foreach f,$(SRCS), mint -q $(MACROS) $(f);)
	@(echo "Done!";)

testall: $(ALLTESTARGS)

sanitytest: $(SANITYTARGS)

##	@(cd MultivariatePowerSeries/PowerSeries/tst ; maple2022 -A 2 -B -b MultivariatePowerSeres.mla fromPolynomial.tst)

mla: 
	rm -f $(MLAFILE);
	cd ..; maple -q MultivariatePowerSeries/build_me.mpl;
	echo "Done!";

clean:
	rm -f $(MLAFILE)