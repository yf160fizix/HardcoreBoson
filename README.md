
# Hardcore Boson Model

## Project Description

Monte Carlo method for the hard core boson model via worm algorithm.

## Model Description
Consider bosons hopping on a periodic cubic lattice whose sites are labeled as $\mathbf{n} = n_x \mathbf{i} + n_y \mathbf{j} + n_z \mathbf{k}$, where $n_x, n_y, n_z = 0, 1, 2, ..., N-1$. The Hamilton operator in Fock space is given by  

$$ H\left( \mu,U \right) = -t \sum\limits_{\mathbf{n},\hat{\alpha}}\left(a_{\mathbf{n}}^{\dagger} a_{\mathbf{n} + \hat{\alpha}} + a_{\mathbf{n}+ \hat{\alpha}}^{\dagger} a_{\mathbf{n} }\right) + (6t -\mu) \sum\limits_{\mathbf{n}} a_{\mathbf{n}}^{\dagger} a_{\mathbf{n}} + U \sum\limits_{\mathbf{n}} a_{\mathbf{n}}^{\dagger} a_{\mathbf{n}}\left( a_{\mathbf{n}}^{\dagger} a_{\mathbf{n}} -1   \right) $$ 

Here we will assume $U$ approaches to infinity so that the bosons will be described as hardcore bosons.

## Executing program

* How to run the program:
```
./run.sh
```
