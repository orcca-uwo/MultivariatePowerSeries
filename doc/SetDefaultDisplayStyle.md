##ENCODING ISO-8859-1
##PROCEDURE(help,nospec,label="MultivariatePowerSeries/SetDefaultDisplayStyle",versionnew="{XENOMORPH}") MultivariatePowerSeries[SetDefaultDisplayStyle]
##TITLE(halfline="Set global display style for power series and univariate polynomial over power series")
##    MultivariatePowerSeries[SetDefaultDisplayStyle]
##ALIAS SetDefaultDisplayStyle, MultivariatePowerSeries:-SetDefaultDisplayStyle, MultivariatePowerSeries
##AUTHOR Ali Asadi masadi4@uwo.ca, Alex Brandt abrandt5@uwo.ca, Marc Moreno Maza moreno@csd.uwo.ca
##
##CALLINGSEQUENCE
##- SetDefaultDisplayStyle('d')
##
##PARAMETERS
##- 'd' : list of equation(s) of the form 'option' = 'value' 
##  where 'option' can be 'maxterms', 'precision', or 'maxdegree', 
##  and 'value' is of type either 'nonnegint' or 'infinity'
##
##DESCRIPTION
##- ~SetDefaultDisplayStyle(d)~ sets the default display style which is used by "Display".
##- The option ~maxterms = n~ sets the maximum number of
##  terms  to be displayed while the the option ~precision = n~
##  sets the maximum degree of the displayed terms.
##- The option 'maxdegree' is used only when displaying a univariate polynomial over power series 'u'. Setting ~maxdegree = n~
##  sets the maximum degree  of a
##  displayed term with respect to the main variable of 'u' only.
##- Using 'infinity' for 'n' indicates that no limits
##  for the number of terms, or their degree, is set.
##- If the default display style is not configured then the display style options 
## will be 'maxterms' = '50', 'precision' = 'infinity', and 'maxdegree' = 'infinity'.
##- The default display style can be overridden for individual power series or univariate polynomials over power series with the "SetDisplayStyle" command. In turn, these settings can be overridden by using options for the "Display" command.
##
##EXAMPLES
##> with(MultivariatePowerSeries):
##> a := GeometricSeries([x, y]);
##> UpdatePrecision(a, 10);
##> b := SumOfAllMonomials([x]);
##> UpdatePrecision(b, 10);
##> SetDefaultDisplayStyle([precision=5]);
##> a; 
##> Display(b);
##> SetDefaultDisplayStyle([maxterms=100, precision=infinity, maxdegree=infinity]);
##> a; 
##
##SEEALSO
##- "Display"
##- "SetDisplayStyle"
##- "GeometricSeries"
##- "SumOfAllMonomials"
##- "HomogeneousPart"
##- "MainVariable"
##
##XREFMAP
##- Display : Help:MultivariatePowerSeries[Display]
##- SetDisplayStyle : Help:MultivariatePowerSeries[SetDisplayStyle]
##- GeometricSeries : Help:MultivariatePowerSeries[GeometricSeries]
##- SumOfAllMonomials : Help:MultivariatePowerSeries[SumOfAllMonomials]
##- HomogeneousPart : Help:MultivariatePowerSeries[HomogeneousPart]
##- MainVariable : Help:MultivariatePowerSeries[MainVariable]
