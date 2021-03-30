##ENCODING ISO-8859-1
##PROCEDURE(help,nospec,label="MultivariatePowerSeries/MultivariatePowerSeries") MultivariatePowerSeries
##TITLE Overview of the MultivariatePowerSeries Package
##ALIAS PowerSeries, MultivariatePowerSeries:-PowerSeries, MultivariatePowerSeries:-UnivariatePolynomialOverPowerSeries
##AUTHOR Ali Asadi masadi4@uwo.ca, Alex Brandt abrandt5@uwo.ca, Marc Moreno Maza moreno@csd.uwo.ca
##CALLINGSEQUENCE
##- MultivariatePowerSeries['command']('arguments')
##- 'command'('arguments')
##
##DEFINE DOTDOTDOT($name,$last)
##INLINEMATH(maple,nospacebefore,notrailingspace,alt="$name[1], ..., $name[$last]")
##  $name[1], `...`, $name[$last]
##ENDINLINEMATH
##ENDDEFINE
##
##DESCRIPTION
##- The `MultivariatePowerSeries` package is a  collection of
##  commands for manipulating multivariate power series and univariate
##  polynomials over multivariate power series.
##- The main algebraic functionalities of this package deal with
##  arithmetic operations (addition, multiplication, inversion, evaluation),
##  for both multivariate power series and univariate
##  polynomials over multivariate power series, as well as
##  factorization of such polynomials.
##
##- Every power series 'q' is encoded by an object storing the following information. 
##  First, a procedure, called the power series *generator*, which, given a non-negative integer
##  'd' returns all non-zero terms of 'q' with total degree 'd'. 
##  Second, a non-negative integer, called the *precision* of 'q', which ensures that all non-zero terms 
##  of 'q' of degree less or equal to that precision have been computed and are stored. 
##  Third, an Array, called the *data array* of 'q'  such that, if
##  all non-zero terms of 'q' of degree 'i' have been computed
##  then they are stored at position 'i' of that array, for all non-negative integers 'i'.
##
##- The implementation of every arithmetic operation, such as addition, multiplication, inversion
##  builds the resulting power series (sum, product or inverse) by creating its generator from the
##  generators of the  operands, which are called *ancestors* of the resulting power series.
##  The coefficients of that resulting power series are computed
##  only when truly needed. Once computed, they are stored in the data array of that power series,
##  where they can be retrieved next time needed.
##  When more terms (than those already stored in the data array) are needed,
##  then the generator is
##  invoked which, in turn, may invoke the generators  of the ancestors.
##
##- The implementation of the factorization commands "WeierstrassPreparation"
##  and "HenselFactorize" is also based on lazy evaluation (also known as calls-by-need).
##  Each factorization command
##  returns the factors as soon as enough information is discovered
##  for initializing the data structures of the factors. The precision
##  of each returned factor, that is, the common precision of its coefficients
##  (which are power series) is  zero. However the generator of each
##  coefficient is known and, thus, the computation of more coefficients
##  can be resumed when a higher precision is requested.
##
##- The commands "PowerSeries" and "UnivariatePolynomialOverPowerSeries" 
##  create power series and univariate
##  polynomials over multivariate power series from objects like polynomials
##  and sequences (given as functions).
##  The commands "GeometricSeries" and "SumOfAllMonomials" create examples
##  of power series.
##- The commands "Display", "SetDefaultDisplayStyle" and "SetDisplayStyle" 
##  control the output format of multivariate power series and univariate
##  polynomials over multivariate power series.
##- The commands "HomogeneousPart", "Truncate", "GetCoefficient", "Precision", "Degree", "MainVariable"
##  access data from a
##  power series or a univariate polynomial over power series.
##- The commands "UpdatePrecision" and "Copy" manipulate data of
##  multivariate power series and univariate
##  polynomials over multivariate power series. 
##- The commands  "Add", "Negate", "Multiply", "Exponentiate", "Inverse", "Divide", "EvaluateAtOrigin"
##  perform arithmetic operations on multivariate power series and univariate
##  polynomials over multivariate power series. The functionality of the first six commands can also be accessed using
##  the standard arithmetic operators when the arguments are power series.
##- The commands "WeierstrassPreparation" and "HenselFactorize" factorize
##  univariate polynomials over multivariate power series. 
##
##INCLUDE assignment_warning.mi
##
##REFERENCES
##-(lead=indent) Alexander Brandt, Mahsa Kazemi, Marc Moreno Maza
##  \"Power Series Arithmetic with the BPAS Library.\"
##  **Computer Algebra in Scientific Computing (CASC)**, **Lecture Notes in Computer Science - 12291**, (2020): 108-128.
##-(lead=indent) Mohammadali Asadi, Alexander Brandt, Mahsa Kazemi, Marc Moreno Maza, and Erik Postma: 
## \" Multivariate Power Series in Maple.\" **Maple Conference 2020, Waterloo, Ontario, Canada, November 2-6, 2020**,
## **Communications in Computer and Information Science (CCIS) series - Springer 2020** (submitted).
##
## 
##SEEALSO
##- "PowerSeries"
##- "UnivariatePolynomialOverPowerSeries"
##- "Display"
##- "SetDefaultDisplayStyle"
##- "SetDisplayStyle"
##- "GeometricSeries"
##- "SumOfAllMonomials"
##- "HomogeneousPart"
##- "Precision"
##- "UpdatePrecision"
##- "GetCoefficient"
##- "GetAnalyticExpression"
##- "Degree"
##- "Copy"
##- "IsUnit"
##- "ApproximatelyEqual"
##- "ApproximatelyZero"
##- "Add"
##- "Negate"
##- "Multiply"
##- "Exponentiate"
##- "Inverse"
##- "Divide"
##- "MainVariable"
##- "Truncate"
##- "EvaluateAtOrigin"
##- "WeierstrassPreparation"
##- "HenselFactorize"

##INDEXPAGE index[package], MultivariatePowerSeries, commands for manipulating multivariate power series and univariate polynomials over multivariate power series

##XREFMAP

##-  PowerSeries : Help:MultivariatePowerSeries[PowerSeries]
##-  UnivariatePolynomialOverPowerSeries : Help:MultivariatePowerSeries[UnivariatePolynomialOverPowerSeries]
##-  Display : Help:MultivariatePowerSeries[Display]
##-  SetDefaultDisplayStyle : Help:MultivariatePowerSeries[SetDefaultDisplayStyle]
##-  SetDisplayStyle : Help:MultivariatePowerSeries[SetDisplayStyle]
##-  GeometricSeries : Help:MultivariatePowerSeries[GeometricSeries]
##-  SumOfAllMonomials : Help:MultivariatePowerSeries[SumOfAllMonomials]
##-  HomogeneousPart : Help:MultivariatePowerSeries[HomogeneousPart]
##-  Precision : Help:MultivariatePowerSeries[Precision]
##-  UpdatePrecision : Help:MultivariatePowerSeries[UpdatePrecision]
##-  GetCoefficient : Help:MultivariatePowerSeries[GetCoefficient]
##-  GetAnalyticExpression : Help:MultivariatePowerSeries[GetAnalyticExpression]
##-  Degree : Help:MultivariatePowerSeries[Degree]
##-  Copy : Help:MultivariatePowerSeries[Copy]
##-  IsUnit : Help:MultivariatePowerSeries[IsUnit]
##-  ApproximatelyEqual : Help:MultivariatePowerSeries[ApproximatelyEqual]
##-  ApproximatelyZero : Help:MultivariatePowerSeries[ApproximatelyZero]
##-  Add : Help:MultivariatePowerSeries[Add]
##-  Negate : Help:MultivariatePowerSeries[Negate]
##-  Multiply : Help:MultivariatePowerSeries[Multiply]
##-  Exponentiate : Help:MultivariatePowerSeries[Exponentiate]
##-  Inverse : Help:MultivariatePowerSeries[Inverse]
##-  Divide : Help:MultivariatePowerSeries[Divide]
##-  MainVariable : Help:MultivariatePowerSeries[MainVariable]
##-  Truncate : Help:MultivariatePowerSeries[Truncate]
##-  EvaluateAtOrigin : Help:MultivariatePowerSeries[EvaluateAtOrigin]
##-  WeierstrassPreparation : Help:MultivariatePowerSeries[WeierstrassPreparation]
##-  HenselFactorize : Help:MultivariatePowerSeries[HenselFactorize
