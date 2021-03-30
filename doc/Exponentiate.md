##ENCODING ISO-8859-1
##PROCEDURE(help,nospec,label="MultivariatePowerSeries/Exponentiate",versionnew="{XENOMORPH}") MultivariatePowerSeries[Exponentiate]
##TITLE(halfline="Exponentiate a power series or a univariate polynomial over power series")
##    MultivariatePowerSeries[Exponentiate]
##ALIAS Exponentiate, MultivariatePowerSeries:-Exponentiate, MultivariatePowerSeries
##AUTHOR Ali Asadi masadi4@uwo.ca, Alex Brandt abrandt5@uwo.ca, Marc Moreno Maza moreno@csd.uwo.ca
##
##CALLINGSEQUENCE
##- 'p'^'e'
##- 'u'^'n'
##- Exponentiate('p', 'e')
##- Exponentiate('u', 'n')
##
##PARAMETERS
##- 'p' : power series generated by this package
##- 'e' : integer
##- 'u' : univariate polynomial over power series generated by this package
##- 'n' : nonnegative integer
##
##DESCRIPTION
##- The commands ~p^e~ and ~Exponentiate(p,e)~ exponentiate the power series 'p'
##  by raising it to the power 'e'.
##- The commands ~u^n~ and ~Exponentiate(u,n)~ exponentiate the univariate polynomial over power
##  series 'u' by raising it to the power 'n'.
##- Note that power series can be raised to any integer power, whereas univariate polynomials over
##  power series can only be raised to nonnegative integer powers.
##
##INCLUDE assignment_warning.mi
##
##EXAMPLES
##> with(MultivariatePowerSeries):
##- We define a power series, _a_.
##> a := GeometricSeries([x,y]):
##<(verification="type") object
##- We can define _a^4_ in three different ways: using multiplication, using the
##  exponentiation operator, or using the `Exponentiate` command.
##> b := Multiply(a,a,a,a):
##<(verification="type") object
##> c := a^4;
##<(verification="type") object
##> d := Exponentiate(a, 4);
##<(verification="type") object
##- We verify that the homogeneous components of _b_, _c_, and _d_ of degree at most 10 are the
##  same.
##> ApproximatelyEqual(b, c, 10);
##< true
##> ApproximatelyEqual(b, d, 10);
##< true
##- We define a univariate polynomial over power series, _f_.
##> f := UnivariatePolynomialOverPowerSeries((z - 1)*(z - 2)*(z - 3) + x*(z^2 + z), z):
##<(verification="type") object
##- Again, we can define _f^3_ in three different ways. We verify that they give the same result (at
##  least for degrees at most 10).
##> g := f * f * f;
##<(verification="type") object
##> h := Exponentiate(f, 3);
##<(verification="type") object
##> k := f^3;
##<(verification="type") object
##> ApproximatelyEqual(g, h, 10);
##< true
##> ApproximatelyEqual(g, k, 10);
##< true
## 
##SEEALSO
##- "GeometricSeries"
##- "Multiply"
##- "UnivariatePolynomialOverPowerSeries"
##- "ApproximatelyEqual"
## 
##XREFMAP
##- GeometricSeries : Help:MultivariatePowerSeries[GeometricSeries]
##- Multiply : Help:MultivariatePowerSeries[Multiply]
##- UnivariatePolynomialOverPowerSeries : Help:MultivariatePowerSeries[UnivariatePolynomialOverPowerSeries]
##- ApproximatelyEqual : Help:MultivariatePowerSeries[ApproximatelyEqual]