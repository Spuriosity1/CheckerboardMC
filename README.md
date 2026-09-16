# Checkerboard Lattice classical Monte Carlo

> A high performance classical Monte Carlo engine written in pure HolyC.

This is a ground-up implementation of Metropois Monte Carlo on the Checkerboard lattice, using pure pointer arithmetic for greater speed and a custom RNG library including some hand-written x86_64 assembly.

# Installing

First, install a compliant HolyC compiler somewhere in PATH -- get it here: https://holyc-lang.com/install.

```bash
git clone https://github.com/Spuriosity1/CheckerboardMC && cd CheckerboardMC
make
```

This will build two binaries -- `ising` and `random_test`. The latter is a test of the bundled splitmix64 PRNG.



## `ising`
Call with `ising <L> <Thot> <Tcold> <seed> <n_step=100> <n_sweep=100> <n_sample=30> <J=1>`.

Sets up a Checkerboard lattice containing $(2L)^2$ sites, i.e. $L^2$ cubic unit cells. Sweeps from Thot to Tcold in n_step steps, sampling an additional n_sample times at the end of each equilibration. $J$ is the interaction strength, in the "signless" convention

$$H = \sum J \sigma_i \sigma_j$$, for $\sigma_i = ±1$. 

i.e. $J=1$ is antiferromagnetic.

## `random_test`
Called with `random_test <N>` to run N Bernoulli trials. Reports deviation from analytical answer for uniform 0-1 distribution.

# to-do

- [ ] Check that the splitmix64 PRNG works more thoroughly
- [ ] Check that the fancier Xoshiro128 PRNG actually works
- [ ] Reorganise the repo
