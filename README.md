Description
--

This repository contains a mix of Magma, Pari and SageMath code for calculating the endomorphism algebras and rings of Jacobian varieties of curves over number fields and finite fields.

Prerequisites
--
An installation of both Pari, Magma and SageMath, so that all of these are available on the command line, is required to run all of the code. For optimal results, set your Pari stack size to a decent size by for example adding
```
parisize = "4096M"
```
to your `~/.gprc` file. This is an optional improvement.

You should also install
```
https://github.com/pascalmolin/hcperiods
```
and include the path to its spec file in your `.magmarc` file. This is not an optimal improvement; the package will not run without it. At any rate Molin--Neurohr's code makes everything far more stable and reliable, so you will want to install it.

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

You need to include a line in your `~/.sage/init.sage` file like
```
magma.attach_spec('[path-to-spec]/spec')
```
similar to the one alluded to in the Magma installation above; if not, the copy of Magma run by SageMath will not see the names in the package by Molin--Neurohr.

Usage
--

A lot of examples, both in SageMath and purely in Magma, are given in the directory `examples/`. The creation of database files, as well as interaction with the LMFDB, is described in the directory `database/`.

More detailed information
--

A description of the data structures used can be found in the files Dicts.md and Descs.md in the directory `documentation/`.
