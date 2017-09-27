CL2QCD
======

CL2QCD is a Lattice QCD application based on OpenCL, applicable to CPUs and GPUs.
It provides the possibility of producing gauge configurations using different algorithms as well as measuring observables on given configurations.

|     Executable     | Description |
| :----------------: | :---------: |
| `su3heatbath`      | Generation of gauge field configurations for SU(3) Pure Gauge Theory                  |
| `hmc`              | Generation of gauge field configurations for Nf=2 (Twisted Mass) Wilson type fermions |
| `rhmc`             | Generation of gauge field configurations for staggered type fermions                  |
| `inverter`         | Measurements of fermionic observables on given gauge field configurations             |
| `gaugeobservables` | Measurements of gauge observables on given gauge field configurations                 |

CL2QCD has been heavily optimized for AMD GPUs, providing world-class performance, but can also be used on NVIDIA GPUs and x86 CPUs.


:warning: Note for users :bangbang:
-----------------------------------

Unfortunately, due to a migration from a different repository hosting service, the history of CL2QCD had to be fully rewritten.
Please consider to clone it again, if you did before the **25.09.2017**.
We apologise for any inconvenience this may cause.


Installation
------------

See [INSTALL](INSTALL) for installation instructions.


Publications
------------

There have been a range of publications on CL2QCD.
If you use CL2QCD for your work please be so kind to cite one of our papers.

* M. Bach, O. Philipsen, C. Pinke, and A. Sciarra [*&laquo;CL2QCD - Lattice QCD based on OpenCL&raquo;*](http://arxiv.org/abs/1411.5219), in The XXXII International Symposium on Lattice Field Theory, Lattice2014, 2014.
* M. Bach, V. Lindenstruth, O. Philipsen, and C. Pinke, [*&laquo;Lattice QCD based on OpenCL&raquo;*](http://arxiv.org/abs/1209.5942), Comput. Phys. Commun., vol. 184, no. 9, p. 19, Mar. 2013.
* M. Bach, V. Lindenstruth, C. Pinke, and O. Philipsen, [*&laquo;Twisted-Mass Lattice QCD using OpenCL&raquo;*](https://inspirehep.net/record/1297645), in The XXXI International Symposium on Lattice Field Theory, Lattice2013, 2013.
* C. Pinke, O. Philipsen, C. Schäfer, L. Zeidlewicz, and M. Bach, [*&laquo;LatticeQCD using OpenCL&raquo;*](http://arxiv.org/abs/1112.5280), 2011.

Many technical details on CL2QCD and the optimizations included can be found in the PhD Thesis of Matthias Bach.

* M. Bach, [*&laquo;Energy- and Cost-Efficient Lattice-QCD Computations using Graphics Processing Units&raquo;*](http://publikationen.ub.uni-frankfurt.de/frontdoor/index/index/docId/37074), PhD Thesis, Goethe Universistät Frankfurt am Main, 2014.

A more recent overview of the code, with a detailed description about the implementation of the staggered formulation, can be found in &sect;3.2 of the PhD Thesis of Alessandro Sciarra.

* A. Sciarra, [*&laquo;The QCD phase diagram at purely imaginary chemical potential from the lattice&raquo;*](https://github.com/AxelKrypton/PhD_Thesis), PhD Thesis, Goethe Universistät Frankfurt am Main, 2016.


License
-------

CL2QCD itself is licensed under the [GPLv3](LICENSE).


However, it includes some bundled third-party code licensed under different licenses.

* Ranlux is distributed under the terms of the GNU General Public (GPL), but it originally was distributed under the [GPLv2](ranlux/COPYING).
* Ranluxcl is distributed under the terms of an MIT license. For details look into the respective [source file](ranluxcl/ranluxcl.cl).
* Einhard is distributed under the the terms of the GNU General Public (GPL)
* Klepsydra is distributed under the the terms of the GNU General Public (GPL)
* The CPP bindings for OpenCL are copyrighted by Khronos and may be [distributed freely](hacks/cl_hpp/CL/cl.hpp).
