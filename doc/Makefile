SRCS = Add ApproximatelyEqual ApproximatelyZero Copy Degree Display Divide EvaluateAtOrigin Exponentiate GeometricSeries GetAnalyticExpression GetCoefficient HenselFactorize HomogeneousPart Inverse IsUnit MainVariable Multiply Negate PowerSeries Precision SetDefaultDisplayStyle SetDisplayStyle Subtract SumOfAllMonomials TaylorShift Truncate UnivariatePolynomialOverPowerSeries UpdatePrecision WeierstrassPreparation MultivariatePowerSeries

TSTFILES = MultivariatePowerSeries:-Add.tst  MultivariatePowerSeries:-ApproximatelyEqual.tst  MultivariatePowerSeries:-ApproximatelyZero.tst  MultivariatePowerSeries:-Copy.tst  MultivariatePowerSeries:-Degree.tst  MultivariatePowerSeries:-Display.tst  MultivariatePowerSeries:-Divide.tst  MultivariatePowerSeries:-EvaluateAtOrigin.tst  MultivariatePowerSeries:-Exponentiate.tst  MultivariatePowerSeries:-GeometricSeries.tst  MultivariatePowerSeries:-GetAnalyticExpression.tst  MultivariatePowerSeries:-GetCoefficient.tst  MultivariatePowerSeries:-HenselFactorize.tst  MultivariatePowerSeries:-HomogeneousPart.tst  MultivariatePowerSeries:-Inverse.tst  MultivariatePowerSeries:-IsUnit.tst  MultivariatePowerSeries:-MainVariable.tst  MultivariatePowerSeries:-Multiply.tst  MultivariatePowerSeries:-Negate.tst  MultivariatePowerSeries:-PowerSeries.tst  MultivariatePowerSeries:-Precision.tst  MultivariatePowerSeries:-SetDefaultDisplayStyle.tst  MultivariatePowerSeries:-SetDisplayStyle.tst  MultivariatePowerSeries:-Subtract.tst  MultivariatePowerSeries:-SumOfAllMonomials.tst  MultivariatePowerSeries:-TaylorShift.tst  MultivariatePowerSeries:-Truncate.tst  MultivariatePowerSeries:-UnivariatePolynomialOverPowerSeries.tst  MultivariatePowerSeries:-UpdatePrecision.tst  MultivariatePowerSeries:-WeierstrassPreparation.tst
RCHELPROOT = $(PWD)/..

MPLTARGS   = $(SRCS:%=%.mpl)

all: clean mpl test

backup:
	$(foreach var, $(SRCS), cp $(var).md $(MAPLELIB)/../help/RegularChains/`dirname $(var)`/;)

mpl: $(MPLTARGS)

test: 
	(cd $(RCHELPROOT)/mw+tst ; $(foreach var, $(TSTFILES), maple <  $(var);))
clean:
	rm -f $(RCHELPROOT)/mw+tst/*.mw $(RCHELPROOT)/mw+tst/*.tst *~ */*~ *.mw  *.tst

%.mpl:
	(cd `dirname $*` ; mpldoc `basename $*`.md ; mv *tst *mw $(RCHELPROOT)/mw+tst 2>/dev/null || : )

m3push:
	(rsync -ru . moreno@192.168.5.20:/home/alex/src/mapleP4/help/MultivariatePowerSeries)

m3pull:
	(cd .. ; rsync -ru moreno@192.168.5.20:/home/alex/src/mapleP4/help/MultivariatePowerSeries . )
