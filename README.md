Description
--

This repository contains a mix of Magma, Pari and SageMath code for calculating the endomorphism algebras and rings of Jacobian varieties of curves over number fields and finite fields.

PLEASE NOTE: The current version of the repository uses Magma version 2.25, which includes the period algorithms by Molin--Neurohr mentioned below.
If you have an earlier version of Magma, please try checking out commit `7433e5` and following the instructions on including the period functionality for old versions that is described below.

Prerequisites
--
This package depends on [`edgarcosta/MagmaPolred`](https://github.com/edgarcosta/MagmaPolred).
The bundle [`edgarcosta/CHIMP`](https://github.com/edgarcosta/CHIMP) includes all the dependencies and optional packages.

An installation of Magma, Pari and (optionally) SageMath, so that all of these are available on the command line, is required to run all of the code.

Most of the algorithms are written in Magma, whose algebro-geometric and numerical capabilities are essential.
Some of the heavy lifting in the creation of number fields and the recognition of complex numbers as algebraic numbers is outsourced to Pari via [`edgarcosta/MagmaPolred`](https://github.com/edgarcosta/MagmaPolred), whose performance when working with number fields is better than that of Magma.
Finally, SageMath is mainly used for the calculation of Frobenius endomorphisms and as a wrapper, because its Python substrate allows for easier creation and manipulation of data structures.


In order to decompose Jacobians, you will need to similarly install and attach [`JRSijsling/curve_reconstruction`](https://github.com/JRSijsling/curve_reconstruction) and its dependencies.


Old magma versions
---

Older than 2.25
----
If your Magma version predates 2.25, then you should also install [`pascalmolin/hcperiods`](https://github.com/pascalmolin/hcperiods) and include the path to its spec file in your `.magmarc` file, using `AttachSpec` in the same way as in the section on the Magma installation below (but with a different target spec file). This is not an optimal improvement; the package will not run without it. At any rate Molin--Neurohr's code makes numerical integration on curves far more stable and reliable, so you will want to install it.

The same holds for code by Christian Neurohr that enables the computation of period matrices of plane curves. It is available via the dependency [`JRSijsling/RiemannSurfaces`](https://github.com/JRSijsling/RiemannSurfaces), a fork containing small modifications of the magnificent original version at [`christianneurohr/RiemannSurfaces`](https://github.com/christianneurohr/RiemannSurfaces).

Older than 2.23
---
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

We highly recommend this package via the bundle [`edgarcosta/CHIMP`](https://github.com/edgarcosta/).
This will include the right dependencies.

The subdirectory `endomorphisms/magma/` includes code that can be run purely within Magma.
You can enable the functionality of this code in Magma by attaching the `endomorphisms/endomorphisms/magma/spec` file with `AttachSpec`.
To make this independent of the directory in which you find yourself, and to active this on startup by default, you may want to indicate the relative path in your `~/.magmarc` file, by adding the line
```
AttachSpec("~/Programs/endomorphisms/endomorphisms/magma/spec");
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

While delegating commands to Magma, the same prerequisites are needed as those mentioned above. For this reason, please create the file `~/.sage/init.sage` if it does not yet exist or add the following line to it:
```
magma.load('~/.magmarc')
```
This ensures that all relevant packages are loaded when outsourcing to Magma.

Usage
--

Examples, both in Magma and in SageMath, are given in the directory `examples/`. The creation of database files, as well as interaction with the LMFDB, is described in the directory `database/`.

Verbose comments are enabled by
```
SetVerbose("EndoFind", m);
SetVerbose("EndoCheck", n);
```
where `m` and `n` are either `1`, `2`, or `3`. A higher value gives more comments.

More detailed information
--

A further description of the data structures and return values used is given in the directory `documentation/`.

Credits
--

The fast calculation of period matrices of hyperelliptic curves in [`pascalmolin/hperiods`](https://github.com/pascalmolin/hcperiods) is based on:

Pascal Molin and Christian Neurohr  
*Computing period matrices and the Abel-Jacobi map of superelliptic curves*  
Mathematics of Computation, 88 (316) (2017)

The calculation of period matrices of plane quartic curves uses the following work:

Christian Neurohr  
*Efficient integration on Riemann surfaces & applications*  
Ph.D. thesis, Carl-von-Ossietzky-Universität Oldenburg (2018)

Citing this code
--

Please cite the following preprint, as well as the two mentioned above, if this code has been helpful in your research.

Edgar Costa, Nicolas Mascot, Jeroen Sijsling, and John Voight  
*Rigorous computation of the endomorphism ring of a Jacobian*
Mathematics of Computation, 88 (2019), 1303-1339 

DOI: https://doi.org/10.1090/mcom/3373 

Preprint available at [arXiv:1705.09248](https://arxiv.org/abs/1705.09248)
