# MultivariatePowerSeries

The **MultivariatePowerSeries** package, which is among new features of [Maple 2021](https://www.maplesoft.com/products/maple/new_features/), 
is an object-oriented implementation of multivariate formal power series and univariate polynomials over such series (UPoPS). 

This package is written in [Maple](https://www.maplesoft.com/) and provides the ability to create and manipulate multivariate power series
with rational or algebraic number coefficients, as well as UPoPS objects whose coefficients are multivariate power series.
Through lazy evaluation techniques and a careful implementation, this package achieves high performance in comparison 
with software counterparts. These power series and UPoPS are employed in optimized implementations of *Weierstrass Preparation Theorem* 
and factorization of UPoPS via *Hensel’s lemma*. 

This figure shows all methods in Maple:
<p align="center"> 
<img  width="800" src="https://github.com/orcca-uwo/MultivariatePowerSeries/blob/master/figures/MPS.png" alt="MPS package">
</p>

## Installation Guide
This package is shipped with **Maple 2021** and can be reached by `with(MultivariatePowerSeries)` command in Maple. 
However, for those using older versions of Maple, there is another way to install this manually as follows:
1. Clone or Download this repository 
2. Go to the repository folder (i.e. `cd ./MultivariatePowerSeries`)
3. Run `make` in your terminal
4. Add `libname :=  "/package/address", libname:` to `~/.mapleinit`

A successful installation creates a new file `MultivariatePowerSeries.mla` at the same directory. 
Note that `/package/address` is the complete path to where `MultivariatePowerSeries.mla` is located.
To test or remove this package, run `make test` or `make clean`. 

## Package Structure
Both power series and UPoPS types are implemented as Maple classes and named, respectively,
  - PowerSeriesObject ([./PowerSeries/src](https://github.com/orcca-uwo/MultivariatePowerSeries/tree/master/PowerSeries/src))
  - UnivariatePolynomialOverPowerSeriesObject ([./UPoPS/src](https://github.com/orcca-uwo/MultivariatePowerSeries/tree/master/UPoPS/src))
  
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
The Maple worksheet [MPS-demo.mw](https://github.com/orcca-uwo/MultivariatePowerSeries/tree/master/demo) is a demo of this package presenting examples such as,
- creating PowerSeries objects from a polynomial or an anonymous function:
<p align="center"> 
<img width="800" src="https://github.com/orcca-uwo/MultivariatePowerSeries/blob/master/figures/PS.png" alt="Example 1">
</p>

- decomposing a UPoPS object using *Weierstrass Preparation Theorem*:
<p align="center"> 
<img width="800" src="https://github.com/orcca-uwo/MultivariatePowerSeries/blob/master/figures/WPT.png" alt="Example 2">
</p>

- factoring a UPoPS object:
<p align="center"> 
<img width="800" src="https://github.com/orcca-uwo/MultivariatePowerSeries/blob/master/figures/HF.png" alt="Example 3">
</p>

## Credit 
This package is designed and developed by M. Asadi under the supervision of M. Moreno Maza and E. Postma 
follows the lazy evaluation scheme of multivariate power series in the [BPAS](https://github.com/orcca-uwo/BPAS) library 
written in C language and implemented by A. Brandt, M. Kazemi, and M. Moreno Maza. 

To cite this work, please use

```
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
