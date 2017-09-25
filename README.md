CL2QCD
======

CL2QCD is a Lattice QCD application based on OpenCL, applicable to CPUs and GPUs.
It provides the following executables:
 * `hmc` -- Generation of gauge field configurations for Nf=2 (Twisted Mass) Wilson type fermions
 * `su3heatbath` -- Generation of gauge field configurations for SU (3) Pure Gauge Theory
 * `rhmc` -- Generation of gauge field configurations for staggered type fermions
 * `inverter` -- Measurements of fermionic observables on given gauge field configurations
 * `gaugeobservables` -- Measurements of gauge observables on given gauge field configurations

CL2QCD has been heavily optimized for AMD GPUs, providing world-class performance, but can also be used on NVIDIA GPUs and x86 CPUs.

:warning: Note for users :bangbang:
-----------------------------------

Unfortunately, due to a migration from a different repository hosting service, the history of CL2QCD had to be fully rewritten.
Please consider to clone it again, if you did before the **25.09.2017**.
We apologise for any inconvenience this may cause.


Installation
------------

See [INSTALL](INSTALL) for installation instructions.

License
-------

CL2QCD itself is licensed under the GPLv3.
However some bundled third-party code is licensed under different licenses.
See [COPYING](COPYING) for details.

Publications
------------

There have been a range of publications on CL2QCD.
If you use CL2QCD for your work please be so kind to cite one of our papers:

* M. Bach, O. Philipsen, C. Pinke, and A. Sciarra “CL2QCD - Lattice QCD based on OpenCL”, in The XXXII International Symposium on Lattice Field Theory, Lattice2014, 2014, http://arxiv.org/pdf/1411.5219.pdf.
* M. Bach, V. Lindenstruth, O. Philipsen, and C. Pinke, “Lattice QCD based on OpenCL”, Comput. Phys. Commun., vol. 184, no. 9, p. 19, Mar. 2013, http://pos.sissa.it/archive/conferences/187/032/LATTICE%202013_032.pdf.
* M. Bach, V. Lindenstruth, C. Pinke, and O. Philipsen, “Twisted-Mass Lattice QCD using OpenCL”, in The XXXI International Symposium on Lattice Field Theory, Lattice2013, 2013, http://arxiv.org/abs/1209.5942.
* C. Pinke, O. Philipsen, C. Schäfer, L. Zeidlewicz, and M. Bach, “LatticeQCD using OpenCL”, 2011, http://arxiv.org/abs/1112.5280.

Many technical details on CL2QCD and the optimizations included can be found in the PhD Thesis of Matthias Bach:

* M. Bach, “Energy- and Cost-Efficient Lattice-QCD Computations using Graphics Processing Units“, PhD Thesis, Goethe Universistät Frankfurt am Main, 2014, http://publikationen.ub.uni-frankfurt.de/frontdoor/index/index/docId/37074.

