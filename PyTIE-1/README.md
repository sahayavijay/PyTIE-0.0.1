## PyTIE_D
PyTIE_D is a free open source Python package for obtaining all the degree-based topological indices expressions and its numerical values for a multiple layer molecular structure with a single parameter. The package was tested under Windows OS.

## Usage

After successfull installation you can immediately use **PyTIE_D** in Pip.

~~~~~~~~~~~~~~~~~~ {.python .numberLines}
    python
    >>> import PyTIE-D   
    >>> from PyTIE-D import *

            
## Structure

Currently the package consists of one modules: **PyTIE_D**
It contains the class method *topological_indices_Degree* and function *topoexpressd*. Using this function, you can obtain all the degree-based topological indices expressions and its numerical values of multi-layered molecular structures with a single parameter, provided you give it the appropriate inputs.

### Calculates topological indices:

* First Zagreb index
* Second Zagreb index
* Hyper Zagreb index
* Third Zagreb index
* Reduced Zagreb index
* Second modified Zagreb index
* Randic index
* Reciprocal Randic index
* Reduced reciprocal Randic index
* General Randic index
* Atom bond connectivity index
* Geometric Arithematic index
* Harmonic index
* Sum-connectivity index
* Inverse sum index
* Alberston index
* Symmetric division index
* Forgotten index
* Sombor index
* bi-Zagreb index
* Tri-Zagreb index
* Geometric Harmonic index
* Geometric Bi-Zagreb index
* Geometric Tri-Zagreb index
* Harmonic-Geometric index

###As Python module
For any UNIX-like system the installation process is trivial:

    pip install pytie-d 

####Depends on:
math
numpy
sympy

##About

PyTIE package written by Sahaya Vijay J, PhD student of Vellore Institute of Technology (VIT University) under supervision of Dr. S. Roy.