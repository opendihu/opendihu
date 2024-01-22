> [!IMPORTANT]  
> This is the new repository for OpenDiHu!
> OpenDiHu was previously hosted in another [repository](https://github.com/maierbn/opendihu), which is no longer maintained. There you can find OpenDiHu Version 1.4 and previous releases. 

> [!WARNING]
> The CI pipeline for this repository is still under construction. 

# Overview
OpenDiHu is a software framework to solve 1D, 2D, and 3D multi-physics problems in parallel with the Finite Element Method.
It is used in the domain of skeletal muscle simulations: Electrophysiology, contraction, neuro-chemo-electro-mechanics.
Design goals are usability, performance and extensibility.

The software is developed at [SGS](https://www.ipvs.uni-stuttgart.de/departments/sgs/) and [IANS](https://www.ians.uni-stuttgart.de/institute/) at the [University of Stuttgart](https://www.uni-stuttgart.de/en/index.html).

Link to Documentation: https://opendihu.readthedocs.io/en/latest/

# Installation
Refer to the documentation for detailed [installation instructions](https://opendihu.readthedocs.io/en/latest/user/installation.html).

However, if you usually skip instructions, do the following:
```
git clone https://github.com/opendihu/opendihu.git && cd opendihu && make
```
If there are error messages, have a look at the log file `config.log`.

# Literature / How to cite

Currently, the following literature written by the authors is available.

* (Open-access) Benjamin's PhD thesis. This is the most comprehensive literature about OpenDiHu and contains details of the implemented models as well as some documentation about the software architecture.

	[Maier, B. (2021). Scalable Biophysical Simulations of the Neuromuscular System. ](https://arxiv.org/abs/2107.07104)

* (Open-access) Journal paper about mesh generation and simulations:

	[Maier, B, Schulte, M. (2022) Mesh generation and multi-scale simulation of a contracting muscle–tendon complex, Journal of Computational Science](https://www.sciencedirect.com/science/article/pii/S1877750322000023)

* Article about computational performance, in proceedings of the HLRS:

	[Krämer, A.; Maier, B.; Rau, T.; Huber, F.; Klotz, T.; Ertl, T.; Göddeke, D.; Mehl, M.; Reina, G.; Röhrle, O. (2020): Multi-physics multi-scale HPC simulations of skeletal muscles, High Performance Computing in Science and Engineering ’20, Transactions of the High Performance Computing Center, Stuttgart (HLRS) 2020, ed. by Nagel, W.; Kröner, D.; Resch, M., Springer International Publishing, 2021](https://link.springer.com/chapter/10.1007/978-3-030-80602-6_13)

* First conference paper:

	[Maier, B., Emamy, N., Krämer, A., & Mehl, M. (2019). Highly parallel multi-physics simulation of muscular activation and EMG. In COUPLED VIII: proceedings of the VIII International Conference on Computational Methods for Coupled Problems in Science and Engineering (pp. 610-621). CIMNE](https://upcommons.upc.edu/handle/2117/190149)
