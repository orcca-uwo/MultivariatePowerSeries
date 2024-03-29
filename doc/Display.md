##ENCODING ISO-8859-1
##PROCEDURE(help,nospec,label="MultivariatePowerSeries/Display",versionnew="{XENOMORPH}") MultivariatePowerSeries[Display]
##TITLE(halfline="nicely display a power series or univariate polynomial over power series")
##    MultivariatePowerSeries[Display]
##ALIAS Display, MultivariatePowerSeries:-Display, MultivariatePowerSeries
##AUTHOR Ali Asadi masadi4@uwo.ca, Alex Brandt abrandt5@uwo.ca, Marc Moreno Maza moreno@csd.uwo.ca
##
##CALLINGSEQUENCE
##- Display('p', 'd')
##- Display('u', 'd')
##
##PARAMETERS
##- 'p' : power series generated by this package
##- 'u' : univariate polynomial over power series generated by this package
##- 'd' : (optional) list of equation(s) of the form 'option' = 'value' where 'option' can be
##  'maxterms', 'precision', or (if the argument is a univariate polynomial over power series)
##  'maxdegree', and 'value' is a non-negative integer or the symbol 'infinity'
##
##DESCRIPTION
##- ~Display(p,d)~ displays the power series 'p' in the style 
##   given by 'd'. The argument 'd' is a list of equations 
## which contains either  'maxterms=n' or 'precision=m', or both. 
## In those equations 'n' and 'm' are of type either 'nonnegint' or 'infinity'.
##- If 'd' is not provided, then the display style of 'p' is taken from
##  an earlier call of the form ~SetDisplayStyle(p,d)~, if any.
##  Otherwise, the display style of 'p' is taken from
## an earlier call of the form ~SetDefaultDisplayStyle(d)~, if any.
##  Finally, if the given display style does not provide values for one (or both) of the settings,
##  Maple uses the default values 'maxterms=50' and 'precision=infinity'.
##- The attribute 'maxterms' sets the maximum number of
##  terms of 'p' to be displayed while the the attribute 'precision'
##  sets the maximum degree of the displayed terms.
##- Using 'infinity' for either 'n' or 'm' indicates that no limits
##  for the number of terms, or their degree, is set.
##- ~Display(u,d)~ displays the coefficients of the
##  univariate polynomial over power series 'u' in in the style
##   given by 'd'.  The argument 'd' can take the same entries as
##   for ~Display(p, d)~, but additionally, the option 'maxdegree=n' can be used:
##  it limits the maximum degree  of a
##  displayed term with respect to the main variable of 'u' only. Like for power series, 'Display'
##  uses 'd' if provided; or otherwise the argument to an earlier call to ~SetDisplayStyle(u, d)~,
##  if any; or otherwise the argument to ~SetDefaultDisplayStyle(d)~, if any; or finally, the same
##  default values 'maxterms=50' and 'precision=infinity' and additionally
##  'maxdegree=infinity'.
##- The limit on the number of terms enforced by the 'maxterms' option is enforced for all terms of
##  the coefficients of 'u' together, the terms of the coefficient power series are not counted
##  independently. Only nonzero terms are counted.
##  
##- The command 'Display' will only ever display coefficients that were computed before;
##  it does not cause computation of further coefficients.
##
##INCLUDE assignment_warning.mi
##
##EXAMPLES
##> with(MultivariatePowerSeries):
##- Define a power series in _x_ and _y_. Initially, only the constant and linear terms are
##  computed, so these are shown.
##> a := GeometricSeries([x, y]);
##<(verification="type") object
##> Display(a);
##- If we compute more terms (in this case, up to homogeneous degree 10), more terms will be displayed.
##> UpdatePrecision(a, 10);
##- The following calling sequence shows up to 20 terms and up to homogeneous degree 5. For this
##  power series, there are a total of 21 terms of homogeneous degree less than or equal to 5, so
##  it shows all but one of them.
##> Display(a,[maxterms = 20, precision = 5]);
##- If we omit the _maxterms_ parameter, its default value is 50, so the following command shows all
##  terms of homogeneous degree less than or equal to 5.
##> Display(a,[precision = 5]);
##- We can set the _maxterms_ and _precision_ parameters to _infinity_ to show all currently
##  computed terms.
##> Display(a,[maxterms = infinity, precision = infinity]);
##- We define a univariate polynomial over power series. Its coefficients have very few terms
##  precomputed, so `Display` doesn't show much.
##> f := UnivariatePolynomialOverPowerSeries([SumOfAllMonomials([x]), GeometricSeries([x, y])], z):
##<(verification="type") object
##> Display(f);
##- If we increase the precision (computing more terms), then `Display` will show more terms.
##> UpdatePrecision(f, 10):
##> Display(f);
##- The _maxdegree_ option cuts the display off after displaying the constant coefficient of the
##  main variable, _z_.
##> Display(f,[maxdegree = 0, precision = 5]);
##- The "SetDisplayStyle" command makes future calls to `Display` use a given set of parameters.
##> SetDisplayStyle(f,[maxterms = 5]);
##> Display(f);
##
##SEEALSO
##- "SetDefaultDisplayStyle"
##- "SetDisplayStyle"
##- "GeometricSeries"
##- "HomogeneousPart"
##- "Truncate"
##- "MainVariable"
##
##XREFMAP
##- SetDefaultDisplayStyle : Help:MultivariatePowerSeries[SetDefaultDisplayStyle]
##- SetDisplayStyle : Help:MultivariatePowerSeries[SetDisplayStyle]
##- GeometricSeries : Help:MultivariatePowerSeries[GeometricSeries]
##- HomogeneousPart : Help:MultivariatePowerSeries[HomogeneousPart]
##- Truncate : Help:MultivariatePowerSeries[Truncate]
##- MainVariable : Help:MultivariatePowerSeries[MainVariable]
