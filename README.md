Description
--

This repository contains a mix of Magma, Pari and SageMath code for calculating the endomorphism algebras and rings of Jacobian varieties of curves over number fields and finite fields.

Prerequisites
--
An installation of Magma, Pari and SageMath, so that all of these are available on the command line, is required to run all of the code. Most of the algorithms are written in Magma, whose algebro-geometric and numerical capabilities are essential. Some of the heavy lifting in the creation of number fields and the recognition of complex numbers as algebraic numbers is outsourced to Pari, whose performance when working with number fields is better than that of Magma. Finally, SageMath is used for the calculation of Frobenius and as a wrapper, because its Python substrate allows for easier creation and manipulation of data structures.

For optimal results, set your Pari stack size to a decent size by for example adding
```
parisize = "4096M"
```
to your `~/.gprc` file. This is an optional improvement.

You should also install [`pascalmolin/hperiods`](https://github.com/pascalmolin/hcperiods) and include the path to its spec file in your `.magmarc` file. This is not an optimal improvement; the package will not run without it. At any rate Molin--Neurohr's code makes numerical integration on curves far more stable and reliable, so you will want to install it.

Finally, in order to decompose Jacobians, you will need [`JRSijsling`](https://github.com/JRSijsling/curve_reconstruction).

Additional prerequisite for older Magma versions
--
In older version of Magma the file `magma/package/Algebra/AlgQuat/interface.m` had the following as line 145:
```
c := [Trace(theta), Norm(theta)];
```
This should be replaced by
```
cpol := MinimalPolynomial(theta);  
assert Degree(cpol) eq 2;  
c := [Coefficient(cpol,1), Coefficient(cpol, 0)];
```

Magma installation 
--

The subdirectory `endomorphisms/magma/` includes code that can be run purely within Magma. You can load all the Magma specific files by attaching the ``endomorphisms/magma/spec`` file with ``AttachSpec``. For example, if you start your session of Magma inside the git directory, you can do this by typing
```
AttachSpec("endomorphisms/magma/spec");
```

SageMath installation
--

To install the package in SageMath, first clone the repository via
```
git clone https://github.com/edgarcosta/endomorphisms.git
```
then go to the newly created directory and type
```
sage -pip install --user --upgrade .
```
After that, a new package called `endomorphisms` will be available for import in SageMath. Once the package is updated on GitHub, pulling the new changes and running the same command will update your installation.

Usage
--

Examples, both in Magma and in SageMath, are given in the directory `examples/`. The creation of database files, as well as interaction with the LMFDB, is described in the directory `database/`.

More detailed information
--

A description of the data structures used in the SageMath wrapper can be found in the files `Dicts.md` and `Descs.md` in the directory `documentation/`.
