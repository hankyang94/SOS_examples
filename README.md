## Preparation
1. Install Matlab
2. Clone this git repo (make sure to pull the submodule spotless as well)
3. cd into spotless in Matlab, and run `spot_install.m` (we use spotless to define polynomial optimization problems)
4. Install [MOSEK](https://www.mosek.com/downloads/) (we use MOSEK as the SDP solver, make sure to get a free academic license)


## Get started
We provide three examples for using **second-order** moment-SOS relaxation to solve nonconvex polynomial optimization problems, and show that the moment-SOS relaxation almost always solves the nonconvex problem to **global optimality**, with certificates. See [here](https://hankyang.seas.harvard.edu/Semidefinite/Moment.html#the-moment-sos-hierarchy) and references therein if you are interested in knowing more.

Note that the Goemans-Williamson MAXCUT relaxation corresponds to a **first-order** moment-SOS relaxation.

### Example 1: binary quadratic program
In this example, we aim to minimize a random quadratic function over a set of binary (+1 and -1) variables. 

Run `example_bqp.m`.

### Example 2: quartic optimization over the unit sphere
In this example, we aim to minimize a random quartic (degree-four) polynomial over the unit sphere (a vector of dimension d constrained to have unit length).

Run `example_q4s.m`

### Example 3: outlier-robust Wahba problem
In this example, we solve a more complicated polynomial optimization. There are two decision variables. 
- q: a dimension-4 vector that is constrained to have unit norm
- theta: a dimension-N vector whose entries are binary variables (i.e., theta(i) = +1 or -1 for every i=1,...,N)

The total number of variables is therefore N+4.

The objective function is a cubic (degree-three) polynomial in q and theta (which has a meaning in computer vision).

Run `example_quasar.m`
