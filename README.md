Description
-----------

This repository contains a mix of Magma, Pari and SageMath code for calculating the endomorphism algebras and rings of Jacobian varieties of curves over number fields and finite fields.

Installation
------------

An installation of both Magma and SageMath is required to run all of the code. To install the package in SageMath, first clone the repository via
```
git clone https://github.com/edgarcosta/endomorphisms.git
```
then go to the newly created directory and type
```
sage -pip install .
```
After that, a new package called `endomorphisms` will be available for import in SageMath.

Magma standalone package
----------------------

The subdirectory `endomorphisms/magma/` includes code that can be run purely within Magma.
You can load all the Magma specific files by attaching the ``endomorphisms/magma/spec`` file with ``AttachSpec``.
For example, if you start your session of Magma inside the git directory, you can do this by typing
 ```
 AttachSpec("endomorphisms/magma/spec");
 ```

Usage
-----

A lot of examples, both in SageMath and purely in Magma, are given in the directory `examples/`. Interaction with databases in the LMFDB is described in the directory `database/`.

A bug fix
---------

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

Optional other packages
-----------------------

For faster calculations of period matrices, you can install the new package
```
https://github.com/pascalmolin/hcperiods
```
In the examples you can then set the flag `molin_neurohr` to `True`. This will make the code run considerably faster, as well as making it more stable.

More information
----------------

A description of the data structures used can be found in the files Dicts.md and Descs.md in the directory `documentation/`.
