# MultivariatePowerSeries

The **MultivariatePowerSeries** package, which is among new features of [Maple 2021](https://www.maplesoft.com/products/maple/new_features/), 
is an object-oriented implementation of multivariate formal power series and univariate polynomials over such series (UPoPS). 

This package is written in [Maple](https://www.maplesoft.com/) and provides the ability to create and manipulate multivariate power series
with rational or algebraic number coefficients, as well as UPoPS objects whose coefficients are multivariate power series.
Through lazy evaluation techniques and a careful implementation, this package achieves high performance in comparison 
with software counterparts. These power series and UPoPS are employed in optimized implementations of *Weierstrass Preparation Theorem* 
and factorization of UPoPS via *Hensel’s lemma*. 

This figure shows all functionalities in Maple:
<p align="center"> 
<img  width="800" src="https://github.com/orcca-uwo/MultivariatePowerSeries/blob/master/figures/MPS.png" alt="MPS package">
</p>

## Package Structure
Our power series and UPoPS types are implemented as Maple classes and named, respectively,
  - PowerSeriesObject
  - UnivariatePolynomialOverPowerSeriesObject
  
The **MultivariatePowerSeries** package contains these two classes along with 
procedures to construct and manipulate objects of those classes. 
In fact, these additional procedures are used to "hide" the object-oriented nature of
the library behind simple procedure calls. The code below shows an overview of 
this package. 

```Ruby 
MultivariatePowerSeries := module()
option package;
  local PowerSeriesObject,
      UnivariatePolynomialOverPowerSeriesObject;

  # create a power series:
  export PowerSeries := proc(...)

  # create a UPoPS:
  export UnivariatePolynomialOverPowerSeries := proc(...)
  
  # Additional procedures to interface these two classes

  module PowerSeriesObject() 
  option object;
    local hpoly :: Array,
          precision :: nonnegint,
          generator :: procedure;
    # other members and methods
  end module;

  module UnivariatePolynomialOverPowerSeriesObject()
  option object;
    local upoly :: Array, 
          vname :: name;
    # other members and methods
  end module;
end module;
```
## Documentation 
A complete documentation is available on [**Maple Online Help**](https://www.maplesoft.com/support/help) and [here](https://github.com/orcca-uwo/MultivariatePowerSeries/tree/master/doc).

## Examples
To create PowerSeries objects from a polynomial or an anonymous function:
<p align="center"> 
<img width="800" src="https://github.com/orcca-uwo/MultivariatePowerSeries/blob/master/figures/PS.png" alt="Example 1">
</p>

To decompose a UPoPS object using *WeierstrassPreparation*:
<p align="center"> 
<img width="800" src="https://github.com/orcca-uwo/MultivariatePowerSeries/blob/master/figures/WPT.png" alt="Example 2">
</p>

To factor a UPoPS object using *HenselFactorize*:
<p align="center"> 
<img width="800" src="https://github.com/orcca-uwo/MultivariatePowerSeries/blob/master/figures/HF.png" alt="Example 3">
</p>

## Credit 
This package has been designed and developed by M. Asadi under the supervision of M. Moreno Maza and E. Postma 
follows the lazy evaluation scheme of multivariate power series in the [BPAS](https://github.com/orcca-uwo/BPAS) library 
written in C language and implemented by A. Brandt, M. Kazemi, and M. Moreno Maza. 

To cite this work, please use

```Latex
@inproceedings{asadi2020MC,
  author = {Asadi, Mohammadali and 
      Brandt, Alexander and 
      Kazemi, Mahsa and 
      Moreno Maza, Marc and 
      Postma, Erik},
  title = {{M}ultivariate {P}ower {S}eries in {M}aple},
  booktitle = {Maple in Mathematics Education and Research},
  series = {Communications in Computer and Information Science},
  publisher = {Springer (accepted)},
  year = {2021}
}
```
