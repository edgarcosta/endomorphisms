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

You should also install [`pascalmolin/hperiods`](https://github.com/pascalmolin/hcperiods) and include the path to its spec file in your `.magmarc` file, using `AttachSpec` in the same way as in the section on the Magma installation below (but with a different target spec file). This is not an optimal improvement; the package will not run without it. At any rate Molin--Neurohr's code makes numerical integration on curves far more stable and reliable, so you will want to install it.

Upcoming installations with Magma will include code by Christian Neurohr that will enable the computation of period matrices of plane curves; you will be able to include these algorithms in `endomorphisms/magma/heuristic/Periods.m` by uncommenting the line marked by the comment `Add Neurohr's code when it becomes available`.

Finally, in order to decompose Jacobians, you will need [`JRSijsling/curve_reconstruction`](https://github.com/JRSijsling/curve_reconstruction). If you do so, do not forget to once again include the path to its spec file in your `.magmarc` file.

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

The subdirectory `endomorphisms/magma/` includes code that can be run purely within Magma.
You can enable the functionality of this code in Magma by attaching the `endomorphisms/endomorphisms/magma/spec` file with `AttachSpec`. To make this independent of the directory in which you find yourself, and to active this on startup by default, you may want to indicate the relative path in your `~/.magmarc` file, by adding the line
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

A description of the data structures used in the SageMath wrapper can be found in the files `Dicts.md` and `Descs.md` in the directory `documentation/`.

Credits
--

The fast calculation of period matrices of hyperelliptic curves in [`pascalmolin/hperiods`](https://github.com/pascalmolin/hcperiods) is based on:

Pascal Molin and Christian Neurohr  
*Computing period matrices and the Abel-Jacobi map of superelliptic curves*  
Mathematics of Computation, 88 (316) (2017)

The (upcoming!) calculation of period matrices of plane quartic curves will use the following work:

Christian Neurohr  
*Efficient integration on Riemann surfaces & applications*  
Ph.D. thesis, Carl-von-Ossietzky-Universität Oldenburg (2018)

Citing this code
--

Please cite the following preprint if this code has been helpful in your research:

Edgar Costa, Nicolas Mascot, Jeroen Sijsling, and John Voight  
*Rigorous computation of the endomorphism ring of a Jacobian*  
Preprint at [arXiv:1705.09248](https://arxiv.org/abs/1705.09248)
