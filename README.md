# CFDEMInterFoamIB
A CFDEM solver to simulate the interaction between two-phase fluid and granualr materials

1. Before installing this solver, ensure that the CFDEM+LIGGGHTS+OpenFOAM are correctly installed.  
**A Chinese tutorial to install the CFDEM+LIGGGHTS+OpenFOAM:  
[如何安装CFDEM+OpenFOAM+LIGGGHTS(COOL三件套)](https://mp.weixin.qq.com/s?__biz=MzI0NzU1NjcyMg==&mid=2247483674&idx=1&sn=ba114d4d4fc7204d755af07217297135&chksm=e9af7e07ded8f71169f83282fc3aa2f9252c6b65eba827ff93b86ef4c084802c97ac033b6629#rd)
2. To install this solver, you just need to run the "remake" file.

## Introduction to tutorial cases
1. single_sphere  
This case is used to show the interaction between particle and two-pase phase in the settling process. Phenomena such as the cavity, splashing, and back-jet.  
2. Multi-sphere
This case is a validation of using the DEM clump. A clump consists of four overlapping sub-spheres is constructed (the particle template is in the location tutorial/multi_sphere/DEM/data/test). Users should learn how to do multi-sphere modeling in LIGGGHTS.

## Animations of some cases:  
### Sphere settling (corresponds to the tutorial case "single_phere"):  
<img src="https://github.com/shenzhihao/CFDEMInterFoamIB/blob/main/animations/settling.gif" width=80% height=80%>  

### Seepage ( not included in tutorial cases):  
<img src="https://github.com/shenzhihao/CFDEMInterFoamIB/blob/main/animations/seepage10mb.gif" width=80% height=80%>  

### Overtopping (not included in tutorial cases):  
<img src="https://github.com/shenzhihao/CFDEMInterFoamIB/blob/main/animations/overtop.gif" width=80% height=80%>  

## Literature about this solver:  
[Shen, Z., Wang, G., Huang, D., & Jin, F. (2022). A resolved CFD-DEM coupling model for modeling two-phase fluids interaction with irregularly shaped particles. Journal of Computational Physics, 448, 110695.](https://www.sciencedirect.com/science/article/pii/S0021999121005908)

## Notice:  
This solver is still on an early release. Questions about installing, compiling, and tutorials cases will be further replied.
