# Sentinel-1 denoising engine
Author: Peter Q Lee. pqjlee (at) uwaterloo (dot) ca
## Installation

You will need to install rust nightly. 

1. Get rustup. Visit [https://rustup.rs] 
2. Install rust nightly. `rustup toolchain install nightly`
3. Install OpenBLAS.
4. HDFlib-dev
5. cd into this directory. Run `python3 setup.py`

## TODO: explain scs requirement.

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




