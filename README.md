# Sentinel-1 noise floor removal engine
Author: Peter Q Lee. pqjlee (at) uwaterloo (dot) ca
## About
This library provides algorithms to remove noise floor intensity patterns that are prevelant in
Sentinel-1 cross-polarized images in EW mode, and to a lesser degree IW mode. This is provided in
four ways\: a command line interface, a Python interface, a Rust interface, and a C interace.

[][]

The algorithms in this library make use of convex optimization to construct noise floors that can be
subtracted to correct the image. Two main algorithms are implemented that make use of the
information provided in the Sentinel-1 xml files.

 1. A linear noise floor removal method that rescales the default noise floor
    that is provided in each Sentinel-1 product. This is the application of the method in
    P. Q. Lee, L. Xu, D. A. Clausi. 2020. Sentinel-1 additive noise removal from cross-polarization
    extra-wide TOPSAR with dynamic least-squares. Remote Sensing of
    Environment. 248. https://doi.org/10.1016/j.rse.2020.111982 . Currently only available in GRD EW mode images.

 2. A non-linear noise floor removal method that computes the noise floor as a power function of
    the antenna pattern. Parameters are estimated with linear programming. This is the application
    of the method in []. Can be applied to both EW and IW GRD mode images.

## Installation

This software has several prerequisites

1. A C compiler
2. The rust compiler/manager: cargo
As these algorithms are written in rust, you will need to install the appropriate tools in order to
compile it for your computer. You will need to install the nightly rust toolchain by doing the followinng.
a. Get rustup. Visit [https://rustup.rs] 
b. Install rust nightly. `rustup toolchain install nightly`

3. OpenBLAS library [].

4. HDF5lib-dev (Development version)[]
5. Splitting conic-solver library [] by:


After these prerequisites have been installed you can install the software by running
`make`

Depending on your system, you may need to supply some additional arguments. 

If SCS isn't found, you can supply the SCS\_LIB=_path_ and SCS\_INC=_path_ environment arguments to
indicate the location of the paths to the splitting conic solver library. On linux you could do this
like:
`SCS_LIB=path/to/scs/out SCS_INC=path/to/scs/include/ make`


## Usage

```python
import denoise_engine
hv, hh, k = denoise_engine.get_dualpol_data(path_to_zip)
```


## Reference functions

`denoise_engine.get_dualpol_data(path_to_zip)`
> Returns: (_crosspol_, _copol_, _k_)
> Where crosspol is the hv or vh image denoised with my method.
> copol is the hh or vv image.
> k is the denoising coefficients produced by my method

`denoise_engine.get_noscale_data(path_to_zip)`
> Returns: (_crosspol_, _copol_, _k_)
> Where _crosspol_ is the hv or vh image denoised with the ESA method.
> _copol_ is the hh or vv image.
> _k_ is an array of 1's.


`denoise_engine.get_raw_data(path_to_zip)`
> Returns (_crosspol_, _copol_, _noisefield_)
> Where _crosspol_ is the unaltered hv or vh from measurement files
> _copol_ is the hh or vh file
> _noisefield_ is the noisefield generated from ESA.




