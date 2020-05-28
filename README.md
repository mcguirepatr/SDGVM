# The SDGVM (est. 1995) #

The Sheffield Dynamic Vegetation Model of Woodward et al 1995 & Woodward & Lomas 2004. Updated with as yet unpublished revisions. 
* Version (versioning system TBD): master used in TRENDY v8, Friedlingstein et al. (2019) ESSD  

Primary authors:
* Ian Woodward
* Mark Lomas

Contributing authors:
* Anthony Walker
* Tristan Quaife
* Ghislain Picard
* Patrick McGuire


### How do I get set up? ###

* Summary of set up - read included files input.doc (up to date, sdescribes hwo to run model) and sdgvm.doc (somewhat out of date, but decribes model science)
* Configuration - run `make clean` then `make` -- requires gfortran
* Dependencies - FORTRAN compiler (gfortran), some R dependencies for post-processing (inc. netcdf libraries if TRENDY output required)
* Database configuration - need met data sets etc, comes separately, not included
* How to run tests - comes with none

### Contribution guidelines ###
Please preserve original model configurations, do not comment out and replace code. 

* Writing tests
* Code review
* Other guidelines - fork repo and make changes on your own branch, pull requests to 'next' branch

### Who do I talk to? ###

* Repo owner or admin - walkerap@ornl.gov (this is fairly unsupported code, will do what I can to help)
* Other community or team contact - Mark Lomas (University of Sheffield)





