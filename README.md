# QCDNUM fork for fragmentation function evolution 

Working copy of [QCDNUM](https://www.nikhef.nl/~h24/qcdnum/) used for evolving fragmentation functions (FF). 

QCDNUM is developed and maintained by M. Botje.

> M. Botje, *Comput. Phys. Commun.* 182 (2011) 490, [arXiv:1005.1481](https://arxiv.org/abs/1005.1481)

Original source and documentation: https://www.nikhef.nl/~h24/qcdnum/



## Fragmentation function evolution

`testjobs/kkD0.f`, `kkD0.cpp`, `kkD0_scan.f` — time-like (fragmentation
function) evolution of the Kniehl-Kramer D0 meson fragmentation function at
LO, following:

> B.A. Kniehl and G. Kramer, "Charmed-Hadron Fragmentation Functions from
> CERN LEP1 Revisited", DESY 06-102, [arXiv:hep-ph/0607306](https://arxiv.org/abs/hep-ph/0607306)

The Kniehl-Kramer paper specifies two separate non-perturbative inputs — a
charm-quark FF (Peterson form) at `mu0 = mc` and a bottom-quark FF (power
law) at `mu0 = mb`. Since DGLAP evolution is linear, each is evolved
separately on the same grid and summed:

```
D_total(x,Q) = D_c-evolved(x,Q) + D_b-evolved(x,Q)   for Q >= mb
D_total(x,Q) = D_c-evolved(x,Q)                      for mc <= Q < mb
```

`kkD0_scan.f` runs the same physics as `kkD0.f` but scans `D_c^D0(x,Q2)`
over `x` at several `Q` values and writes the result to a data file for
plotting.

## Build & run

Full install instructions (Autotools, local build, OpenMP) are in
[`README`](README). Quick start:

```bash
./configure
make                # or make -j8
make install        # or sudo make install

cd run
./runtest kkD0.f        # or kkD0.cpp / kkD0_scan.f
```

## License

QCDNUM is distributed under the GNU General Public License v3 — see
[`COPYING`](COPYING) and [`LICENCE`](LICENCE). Custom additions in this
repo follow the same terms.
