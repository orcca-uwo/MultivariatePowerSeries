##ENCODING ISO-8859-1
##PROCEDURE(help,nospec,label="MultivariatePowerSeries/Inverse",versionnew="{XENOMORPH}") MultivariatePowerSeries[Inverse]
##TITLE(halfline="Compute the inverse of a power series")
##    MultivariatePowerSeries[Inverse]
##ALIAS Inverse, MultivariatePowerSeries:-Inverse, MultivariatePowerSeries
##AUTHOR Ali Asadi masadi4@uwo.ca, Alex Brandt abrandt5@uwo.ca, Marc Moreno Maza moreno@csd.uwo.ca
##
##CALLINGSEQUENCE
##- 1/'p'
##- Inverse('p')
##
##PARAMETERS
##- 'p' : power series generated by this package
##
##DESCRIPTION
##- The commands ~1/p~ and ~Inverse(p)~ compute the multiplicative inverse of the power series
##  'p'. This requires that 'p' is invertible, that is, that 'p' has a non-zero constant term; if
##  that is not the case, an error is signaled.
##
##INCLUDE assignment_warning.mi
##
##EXAMPLES
##> with(MultivariatePowerSeries):
##- We define a power series, _a_, representing a polynomial.
##> a := PowerSeries(1-x-y);
##<(verification="type") object
##- We can define its inverse in two equivalent ways:
##> b := 1/a;
##<(verification="type") object
##> c := Inverse(a);
##<(verification="type") object
##- We verify that the two definitions are equal, at least for the terms up to homogeneous degree
##  10.
##> Truncate(b - c, 10);
##< 0
##- A different power series represents the sine of _x_.
##> sx := PowerSeries(d -> ifelse(d :: odd, (-1)^((d-1)/2) * x^d/d!, 0), analytic = sin(x));
##<(verification="type") object
##- Because the constant coefficient of _sx_ is zero, we cannot invert it. (Its multiplicative
##  inverse is a Laurent series, not a power series.)
##> 1/sx;
##<(verification="testerror") "not invertible"
##
##SEEALSO
##- "PowerSeries"
##- "GeometricSeries"
##- "Divide"
##- "Inverse"
##- "IsUnit"
##- "ApproximatelyEqual"
##
##XREFMAP
##- PowerSeries : Help:MultivariatePowerSeries[PowerSeries] 
##- GeometricSeries : Help:MultivariatePowerSeries[GeometricSeries] 
##- Divide : Help:MultivariatePowerSeries[Divide]
##- Inverse : Help:MultivariatePowerSeries[Inverse]
##- IsUnit : Help:MultivariatePowerSeries[IsUnit]
##- ApproximatelyEqual : Help:MultivariatePowerSeries[ApproximatelyEqual]
