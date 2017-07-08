# Proceed with caution

We are in the process of merging two repositories:

* https://github.com/JRSijsling/heuristic_endomorphisms
* https://github.com/edgarcosta/genus2_endomorphisms

and some functionality might be broken.
# Description

This repository contains a mix of Magma, Pari and Sage code for calculating a heuristic (and usually correct) approximation of the endomorphism algebras and rings of Jacobian varieties of hyperelliptic curves.
For now, the code assumes the curves involved to be defined over QQ, and in many cases assumes that we are in genus 2.
The repository also provides a frame work to rigorously verify the above computations.

# Optional other packages

If you have the (not yet publically available) Magma code for calculating period matrices by Pascal Molin and Christian Neurohr, then please attach it in your `~/.magmarc` file. You can then set the flag `have_oldenburg` to `True`, which will make the code run considerably faster, as well as making it more stable.

# Installation

An installation of both Magma and Sage is required to run this code. 

To install the package pick a directory where you would like to place a copy of the package and do:
```
git clone https://github.com/edgarcosta/endomorphisms.git
```

# Loading the package
The `endomorphisms` package can be loaded in Sage in three ways.
Let `[PATH]` denote the path for the directory into which you have cloned or copied the repository.

### 1. As a python module
```
sys.path.append('[PATH]')
from endomorphisms import *
``` 
to make it more convenient you could also add the line
```
sys.path.append('[PATH]')
```
the Sage initialization file (typically found in  `~/.sage/init.sage`).

### 2. Using `load()`
By going to `[PATH]/endomorphisms` directory and typing
```
load('Initialize.sage')
```

### 3. Loaded at startup 
If you prefer, you can have Sage to load every time you start Sage by adding the following lines to the Sage initialization file (typically found in ~/.sage/init.sage):
```
__endodir__ = '[PATH]/endomorphisms'  
load(__endodir__ + 'Initialize.sage')
```

alternatively

```
sys.path.append('[PATH]')
from heuristic_endomorphisms import *
```

Note that this will also startup Magma automatically.

# Usage 

A sample of examples are given in `examples/`.

# A bug fix

It is highly recommended to fix a Magma bug before using this package. In old version the file `magma/package/Algebra/AlgQuat/interface.m` had the following as line 145:
```
c := [Trace(theta), Norm(theta)];
```
This should be replaced by
```
cpol := MinimalPolynomial(theta);  
assert Degree(cpol) eq 2;  
c := [Coefficient(cpol,1), Coefficient(cpol, 0)];
```

# Data structures used

The most important data returned by the code are the classes `OverField`, `Lattice` and `Decomposition` that can be created from the main class `EndomorphismData`. Given an instance `End` of endomorphism data, these classes are formed by the function evaluations `End.over_field(K)`, `End.lattice()`, and `End.decomposition()`. We briefly describe the corresponding data structures.

In all cases, the first evaluation uses Magma to calculate a raw version of all relevant data. At first, the user is only presented with a raw form of part of this data, consisting of a recursive list of strings and integers. This raw form can be accessed by a call of the `description` method, and a pretty print form is returned when calling `pretty_print`.

This raw form does not contain all objects of interest. For example, it does not give the endomorphism ring as an explicit algebra on generators. The user can ask for more detailed information using various other submethods, for which we refer to `Wrapper.sage`. When asking for such data for the first time, Sage will create a recursive dictionary with all information that the user could ever request. Here is a description of the various dictionaries thus returned:

1. `OverField`: This is a dictionary with three keys, namely `representation`, `algebra`, and `description`. Because of its importance for quicker input-ouput routines, we describe the `description` dictionary in more detail in the next entry. The other two entries are as follows:
* `representations` gives a list of dictionaries, which corresponds to the generators of the endomorphism ring in their various representations. These dictionaries have the keys `tangent`, `homology`, `approx` and `corresp`. The first three of these give various linear representations; `tangent` gives the action on the tangent space of the Jacobian, of which `approx` is a complex approximation, and `homology` the action on the homology. The final entry `corresp` corresponds to a corresponding correspondence. This will only be calculated when needed or requested, since it requires more time.
* `algebra` gives another dictionary with keys `alg_QQ`, `alg_ZZ`,  `alg_RR`, and `alg_ST`. The key `alg_QQ` is a Magma object isomorphic with the endomorphism algebra, and the key `alg_ZZ` gives generators for the endomorphism ring as a subring of this endomorphism algebra. The key `alg_RR` is a list of strings that describes the tensor product of the endomorphism ring with RR, while `alg_ST` gives a string that represents the Sato-Tate group.

2. The description key is the only of these dictionaries whose pendant is an object in Sage (not in Magma). This object is nothing but a recursive list of strings and integers. The corresponding dictionary instead has four keys, namely `factors_QQ`, `desc_ZZ`, `desc_RR` and `sato_tate`. The key `desc_RR`, as `alg_RR` above, is a list of strings that describes the tensor product of the endomorphism ring with RR, while `sato-tate` gives a single string describing the Sato-Tate group (when the classification is known; otherwise it is the string `undef`). The other two keys are once more composed objects:
* `factors_QQ` is a list of dictionaries. The dictionaries in this list correspond to the factors of the endomorphism algebra after central decomposition. They have keys `albert_type`, `base_field`, `dim_sqrt` and `disc`. The key `albert_type` describes the extended Albert classification of the endomorphism algebra by a string. The key `base_field` describes the center of this algebra, as a list of integers. The key `dim_sqrt` is an integer that is the square root of the dimension of the corresponding algebra over this field, and `disc` gives the norm of its discriminant, if applicable.
* `desc_ZZ` is a much smaller dictionary and describes some of the structure of the endomorphism ring proper. Its keys are `index` and `is_eichler`. The first of these gives the index of the endomorphism ring in a maximal order of the endomorphism algebra that contains it, and the second indicates whether or not this order is Eichler.

3. `Lattice`: A list of dictionaries with keys `field` and `structure`. The key `field` has a value that is itself a dictionary with key `seq` and `magma`. The key `seq` describes the number field via a string of integers (defining a corresponding minimal polynomial), and `magma` gives the corresponding Magma object. The key `structure` in the dictionary is a dictionary as used by OverField.

4. `Decomposition`: A list of dictionaries with keys `field`, `idem`, `factor`, and `proj`; these correspond the factors in the decomposition of the Jacobian. The key `field` gives the base field of the corresponding factor and projection; it in principle superfluous but included for clarity. The key `idem` describes the corresponding idempotent endomorphism via the same dictionary used for the elements of `gens` in `OverField`. The key `factor` is also a dictionary with two representations of this factor, namely `analytic` and `algebraic`. Finally, the key `proj` is like `idem`, but describes the projection to the factor instead of the idempotent itself.
