# Dynamic Bandwidth Variational Mode Decomposition

**Authors**: Andreas G. Angelou, Georgios K. Apostolidis, Leontios J. Hadjileontiadis

**Brief introduction**: Dynamic Bandwidth Variational Mode Decomposition (DB-VMD) is a signal decomposition method and a generalization of Variational Mode Decomposition (VMD). In particualr, DB-VMD addresses the constant bandwidth Wiener filters limitation of VMD and proposes a scheme for Wiener filters with dynamic bandwidth. Experiments in synthetic signals underscore DB-VMD’s superior noise robustness and adaptability in comparison to VMD, paving the way for many applications, especially when the analyzed signals are contaminated with noise.

## The repository is structured as follows:
- *Method_Scripts* folder: Script implementations of DB-VMD and VMD
- *Experiments* folder: Experiment scripts about **tone seperation** and **noise robustness**. Also, a script named **visualize_methods.m** is added for visualizing DB-VMD and VMD

### User should run the scripts included in the *Experiments* folder.


**Acknowledgments: The VMD Matlab code was taken from the publicly available MATLAB Central File Exchange page below:**
VMD [1]: Dominique Zosso (2021). Variational Mode Decomposition [MATLAB Central File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/44765-variational-mode-decomposition).

The corresponding paper:
[1] [K. Dragomiretskiy and D. Zosso, "Variational Mode Decomposition," in IEEE Transactions on Signal Processing, vol. 62, no. 3, pp. 531-544, Feb.1, 2014, doi: 10.1109/TSP.2013.2288675](https://ieeexplore.ieee.org/abstract/document/6655981).
