SRCS = MultivariatePowerSeries/src/MultivariatePowerSeries.mpl
MLAFILE = MultivariatePowerSeries.mla
ASSERTLEVEL = 2

ALLTESTFILES=MultivariatePowerSeries/PuiseuxSeries/tst/addition.tst MultivariatePowerSeries/PuiseuxSeries/tst/convert_to_pso.tst MultivariatePowerSeries/PuiseuxSeries/tst/exports.tst MultivariatePowerSeries/PuiseuxSeries/tst/inverse.tst MultivariatePowerSeries/PuiseuxSeries/tst/longerExample1.tst MultivariatePowerSeries/PuiseuxSeries/tst/longerExample2.tst MultivariatePowerSeries/PuiseuxSeries/tst/makeCompatible.tst MultivariatePowerSeries/PuiseuxSeries/tst/multiplication.tst MultivariatePowerSeries/PuiseuxSeries/tst/objectDefinition.tst MultivariatePowerSeries/PuiseuxSeries/tst/randomtest_MakeRaysCompatible.tst MultivariatePowerSeries/PuiseuxSeries/tst/randomtest_objectDefinition.tst MultivariatePowerSeries/PuiseuxSeries/tst/trigonometricIdentity1.tst MultivariatePowerSeries/PuiseuxSeries/tst/trigonometricIdentity2.tst MultivariatePowerSeries/PuiseuxSeries/tst/truncate.tst MultivariatePowerSeries/PowerSeries/tst/addition.tst MultivariatePowerSeries/PowerSeries/tst/constant.tst MultivariatePowerSeries/PowerSeries/tst/deepCopy.tst MultivariatePowerSeries/PowerSeries/tst/FactorOutMonomial.tst MultivariatePowerSeries/PowerSeries/tst/fromPolynomial.tst MultivariatePowerSeries/PowerSeries/tst/fromProcedure.tst MultivariatePowerSeries/PowerSeries/tst/geometricSeries.tst MultivariatePowerSeries/PowerSeries/tst/identity.tst MultivariatePowerSeries/PowerSeries/tst/inverse.tst MultivariatePowerSeries/PowerSeries/tst/naryAdd.tst MultivariatePowerSeries/PowerSeries/tst/naryMul.tst MultivariatePowerSeries/PowerSeries/tst/one.tst MultivariatePowerSeries/PowerSeries/tst/product.tst MultivariatePowerSeries/PowerSeries/tst/Substitute.tst MultivariatePowerSeries/PowerSeries/tst/Substitute-wbv.tst MultivariatePowerSeries/PowerSeries/tst/sumOfAllMonomials.tst MultivariatePowerSeries/PowerSeries/tst/TaylorShift.tst MultivariatePowerSeries/PowerSeries/tst/zero.tst MultivariatePowerSeries/UPoPS/tst/addition.tst MultivariatePowerSeries/UPoPS/tst/constant.tst MultivariatePowerSeries/UPoPS/tst/fromPolynomial.tst MultivariatePowerSeries/UPoPS/tst/hensel_factorization.tst MultivariatePowerSeries/UPoPS/tst/one.tst MultivariatePowerSeries/UPoPS/tst/product.tst MultivariatePowerSeries/UPoPS/tst/Puiseux_theorem.tst MultivariatePowerSeries/UPoPS/tst/taylorShift.tst MultivariatePowerSeries/UPoPS/tst/UPoPSOverPuSO.tst MultivariatePowerSeries/UPoPS/tst/wp.tst MultivariatePowerSeries/UPoPS/tst/zero.tst

SOMETESTFILES = MultivariatePowerSeries/PowerSeries/tst/addition.tst MultivariatePowerSeries/PowerSeries/tst/constant.tst MultivariatePowerSeries/PowerSeries/tst/deepCopy.tst MultivariatePowerSeries/PowerSeries/tst/FactorOutMonomial.tst MultivariatePowerSeries/PowerSeries/tst/fromPolynomial.tst MultivariatePowerSeries/PowerSeries/tst/fromProcedure.tst MultivariatePowerSeries/PowerSeries/tst/geometricSeries.tst MultivariatePowerSeries/PowerSeries/tst/identity.tst MultivariatePowerSeries/PowerSeries/tst/inverse.tst MultivariatePowerSeries/PuiseuxSeries/tst/addition.tst

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