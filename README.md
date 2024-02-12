# Readme for reconstructing curves of genus 4
This repository contains code to reconstruct a genus 4 curve from its theta constants. See https://arxiv.org/abs/2402.03160 for the accompanying article.

## Structure of the repository
The repo contains two folders:

- magma
- examples

In the Magma folder we have the folowing files:
- **auxiliary.m** Contains a few auxiliary methods.
- **azygetic.m** Contains the method that takes a system of four odd azygetic characteristics and extends it to a fundamental system.
- **fast_theta.m** An adaption of Labrande's fast theta code to compute Theta null values.
- **gluingfuncs.m** Contains an array of maximal isotropic subgroups of F_2^4 and a few methods that are helpful for gluing.
- **good-coordinates.m** Contains tools for arithmetic construction from small period matrices
- **igusa_quartic.m** Reconstruct a genus 2 curve from the fourth powers of the thetas.
- **reconstruction.m** The main algorithm to reconstruct a genus 4 curve from its theta constants.
- **rosenhain.m** Compute the Rosenhain invariants of a genus 4 curve from the squares of its theta constants.
- **schottky.m** Computes the value of the Schottky modular form at a given period matrix tau.
- **signs.m** Code that uses Riemann's formula to correct the signs of the thetas that has been lost due to taking squares of thetas.
- **theta.m** Contains methods that compute Theta functions and the Siegel reduction of a period matrix.

In the examples folder we have the following files:
- **noam_example.m** Constructs a genus 4 curve whose ... TODO

## How to run the examples

