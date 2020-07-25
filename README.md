# The SDGVM (est. 1995) #

The Sheffield Dynamic Vegetation Model of Woodward et al 1995 & Woodward & Lomas 2004. Updated with as yet unpublished revisions. 


Primary authors:

* Ian Woodward
* Mark Lomas

Contributing authors:

* Anthony Walker
* Tristan Quaife
* Ghislain Picard
* Patrick McGuire


### How do I get set up? ###

* Summary of set up - read included files in `documentation` directory `input.doc` (up to date, describes how to run model), `projects.docx` (up to date, helps set up ensembles and HPC runs), and `sdgvm.doc` (somewhat out of date, but describes core model science).
* Configuration - run `make clean` then `make` -- requires gfortran
* Dependencies - FORTRAN compiler (gfortran), some R dependencies for post-processing (inc. netcdf libraries if TRENDY output required)
* Database configuration - need met data sets etc, comes separately, not included
* How to run tests - comes with none

### Contribution guidelines ###
Fork repo and make changes on your own branch, make pull requests to `next` branch.
Preserve original model configurations.

* do not comment out and replace code unless a bug, use input file switches for new development
* do not change defaults in code  
* Writing tests - none
* Code review - ad hoc, when I can get to it 
* Other guidelines - 

### Who do I talk to? ###

* Repo owner or admin - walkerap@ornl.gov (this is fairly unsupported code, will do what I can to help)
* Other community or team contact - Mark Lomas (University of Sheffield)





