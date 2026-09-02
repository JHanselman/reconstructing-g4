# Readme for reconstructing curves of genus 4
This repository contains code to reconstruct a genus 4 curve from its theta constants. See https://arxiv.org/abs/2402.03160 for the accompanying article.

## Structure of the repository
The repo contains four folders:

- magma
- examples
- Cext
- tests

In the Magma folder we have the folowing files:
- **auxiliary.m** Contains a few auxiliary methods.
- **azygetic.m** Contains the method that takes a system of four odd azygetic characteristics and extends it to a fundamental system.
- **fast_theta.m** An adaption of Labrande's fast theta code to compute Theta null values.
- **FlintWrapper.m** `ThetaFlint(z, tau : ord)`: theta functions (and their derivatives in z) for all characteristics, computed rigorously by FLINT's `acb_theta` through the external program in `Cext/` (see below).
- **gluingfuncs.m** Contains an array of maximal isotropic subgroups of F_2^4 and a few methods that are helpful for gluing.
- **good-coordinates.m** Contains tools for arithmetic construction from small period matrices
- **igusa_quartic.m** Reconstruct a genus 2 curve from the fourth powers of the thetas.
- **reconstruction.m** The main algorithm to reconstruct a genus 4 curve from its theta constants.
- **rosenhain.m** Compute the Rosenhain invariants of a genus 4 curve from the squares of its theta constants.
- **schottky.m** Computes the value of the Schottky modular form at a given period matrix tau.
- **signs.m** Code that uses Riemann's formula to correct the signs of the thetas that has been lost due to taking squares of thetas.
- **theta.m** Contains methods that compute Theta functions and the Siegel reduction of a period matrix.

The examples folder has two subfolders. In **paper-examples** are the computations from the article:
- **Example-Gluing.m** Glues the Jacobians of two genus 2 curves along their 2-torsion and reconstructs a genus 4 curve over Q with the glued Jacobian.
- **Example-Modular.m** Reconstructs a genus 4 curve over Q from the big period matrix (hard-coded, 300 digits) of the modular abelian fourfold attached to the newform 778.2.a.a; **778.2.a.a.m** is the LMFDB download of that newform, used to verify the curve.
- **QMexample.m**, **RMexample.m** Search the Igusa quartic for genus 2 curves that glue to a genus 4 curve whose Jacobian has quaternionic, resp. real, multiplication. **Galois.m** (Galois action and cocycle on the Igusa coordinates) and **g2-rmlmfdb.m** (genus 2 curves with RM from the LMFDB) are loaded by these.
- **Finite-Field-example.m** Runs the reconstruction from tritangent data over F_37.

In **construct-examples**, **construct-Lemma-Section-3.2-proofs.m** constructs a plane quartic from an Aronhold system of bitangents and **Aronhold-check.m** verifies the identities used in the proofs of Section 3.2.

In the Cext folder, **acb_theta_cli.c** is a small C program linked against FLINT that evaluates theta functions and their derivatives; `ThetaFlint` in magma/FlintWrapper.m calls it. It replaces the earlier python-flint route (CHIMP's Theta.magma) and needs nothing but libflint >= 3.1.

In the tests folder:
- **theta_flint_vs_chimp.m** compares `ThetaFlint` with CHIMP's python-flint values, the pure-Magma theta series, parity of theta constants and finite differences.
- **run_example_gluing.m**, **run_example_modular.m** run the two paper examples that go through the theta code and verify the results independently (see below).

## Building the FLINT helper

```
cd Cext
make          # Linux with libflint-dev >= 3.1
make mac      # macOS: brew install flint mpfr gmp
make test     # self-test
```
See Cext/README.md for a server (no-admin) build. The Magma code finds `Cext/acb_theta_cli` automatically; to use a binary elsewhere set `THETA_FLINT_BIN=/path/to/acb_theta_cli`. If CHIMP is attached as well, note that both packages define `ThetaFlint(z, tau)`; attach this package last so that the FLINT-only version is used.

## How to run the examples

All examples need Magma and the compiled FLINT helper (see the previous section). Attach the package with

```
AttachSpec("magma/spec");
```

adjusting the path as needed; the example files do this themselves with paths relative to `examples/paper-examples`, so run them from that directory.

**Gluing** (needs nothing beyond this repository; about a second):

```
cd examples/paper-examples
magma Example-Gluing.m
```

This computes the period matrices of the two genus 2 curves at 200 digits, the maximal isotropic subgroup `V` given by the identification of their 2-torsion, the glued big period matrix `Q`, and finally `Eqs := RationalReconstructCurveG4(Q)`, a quadric and a cubic in P^3 over Q. To also check the result, run the test driver instead:

```
magma ../../tests/run_example_gluing.m
```

It verifies that the Frobenius traces of the recovered curve agree with those of `Jac(X1) x Jac(X2)` at all good primes up to 200 (up to the sign of a quadratic character, which it identifies: the isogeny is defined over Q(sqrt(-730)); since a non-hyperelliptic curve has no quadratic twists, this is the only curve over Q in the picture) and, if Thomas Bouchet's [Genus-4](https://github.com/Thittho/Genus-4) and [GL-Equivalence](https://github.com/Thittho/GL-Equivalence) packages are found in `~/github` or `~`, that its invariants agree with the model recorded in `QMexample.m`.

**Modular** (needs nothing beyond this repository; a couple of seconds, plus about 20 s for the verification):

```
cd examples/paper-examples
magma ../../tests/run_example_modular.m
```

This reads the 300-digit big period matrix from `Example-Modular.m`, runs `RationalReconstructCurveG4` on it and `ReconstructCurveG4` over C on the Siegel-reduced small period matrix, and then checks, at every good prime p <= 60, that the L-polynomial of the recovered curve is the product over the four embeddings of the Hecke field of `1 - a_p T + p T^2`, using the eigenvalues `a_p` of 778.2.a.a stored in `778.2.a.a.m`. `Example-Modular.m` itself can also be run directly (`magma Example-Modular.m`; it expects `GL-Equivalence` and `Genus-4` to be checked out next to this repository). Its own verification against the newform up to p = 1000, which rebuilds the newform with `ModularForms`, is **switched off by default**: it is guarded by the flag `verify_with_newform := false` near the top of the file, and the example prints a reminder when it is skipped. Set the flag to `true` to run it; expect it to take much longer than the reconstruction.

**QM and RM** (`QMexample.m`, `RMexample.m`) are searches rather than single computations: they loop over rational points on the Igusa quartic and test each candidate gluing, which takes hours. `RMexample.m` also needs the `EndomorphismAlgebra` intrinsic from the [endomorphisms](https://github.com/edgarcosta/endomorphisms) package (included in [CHIMP](https://github.com/edgarcosta/CHIMP)). Run them from `examples/paper-examples` with `magma QMexample.m` or `magma RMexample.m`; `QMexample.m` writes its log to `Output.txt`. The pair of genus 2 curves found by the QM search is the one glued in `Example-Gluing.m`.

**Finite field** and **construct-examples** are self-contained: `magma Finite-Field-example.m` in `examples/paper-examples`, and `magma Aronhold-check.m` in `examples/construct-examples`.

**Testing the theta code** after building or updating FLINT:

```
cd Cext && make test && cd ..
magma tests/theta_flint_vs_chimp.m
```

The second command compares `ThetaFlint` against the pure-Magma theta series, checks the parity and numbering of the characteristics, and checks FLINT's derivatives against finite differences; if a `python3` with `python-flint` is available (or CHIMP's own venv exists), it also compares against the values CHIMP's python-flint wrapper produces.
