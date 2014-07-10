# libmva v0.9

libmva comprises of basic matrix and linear algebra, basic statistics (standard 
deviation, variance, etc.), Principal Component Analysis (PCA) and Projection to
Latent Structures (PLS).  All techniques support handling of missing data where
appropriate.

## Compile
Currently, libmva is compiled with clang and mainly developed on Mac.  It will compile on Linux w/ GCC or clang and Windows with MSVC, but I've not yet created a generic platform independent Makefile. So for now from the root folder (with clang installed) you can simply run:

  make
  
This will build the static lib (libmva.a) in the `mva` folder and the test program (mvlib_test) under the `test` folder.

### Debug

To build with debug symbols (-g) resulting in `libmvad.a` and `mvlib_testd`:
  
  make debug
  
### Clean

Clean all the object files with:

  make clean
  
Clean everything including the build artifacts:

  make distclean

## Third-party licensing notes

libmva uses Cephes library (http://www.netlib.org/cephes/).  This library does 
not contain any specific license. 

# License

libmva is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

