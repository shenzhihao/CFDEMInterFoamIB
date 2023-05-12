# CFDEMInterFoamIB
A CFDEM solver to simulate the interaction between two-phase fluid and granualr materials

1. Before installing this solver, ensure that the CFDEM+LIGGGHTS+OpenFOAM are correctly installed.  
**A Chinese tutorial to install the CFDEM+LIGGGHTS+OpenFOAM:  
[如何安装CFDEM+OpenFOAM+LIGGGHTS(COOL三件套)](https://mp.weixin.qq.com/s?__biz=MzI0NzU1NjcyMg==&mid=2247483674&idx=1&sn=ba114d4d4fc7204d755af07217297135&chksm=e9af7e07ded8f71169f83282fc3aa2f9252c6b65eba827ff93b86ef4c084802c97ac033b6629#rd)
2. To install this solver, you just need to run the "remake" file after you successfully installed the previous codes.

## Introduction to tutorial cases
1. single_sphere  
This case is used to show the interaction between particle and two-pase phase in the settling process. Phenomena such as the cavity, splashing, and back-jet can be observed.  
2. Multi_sphere_fish  
This case is an application of using the DEM clump. A clump consists of overlapping sub-spheres is constructed to represent a fish-shaped object. Users should learn how to do multi-sphere modeling in LIGGGHTS.
<!-- This case is a validation of using the DEM clump. A clump consists of four overlapping sub-spheres is constructed (the particle template is in the location tutorial/multi_sphere/DEM/data/test). Users should learn how to do multi-sphere modeling in LIGGGHTS. -->

## Animations of some cases:  
### Sphere settling (corresponds to the tutorial case "single_sphere"):  
<img src="https://github.com/shenzhihao/CFDEMInterFoamIB/blob/main/animations/settling.gif" width=80% height=80%>  
Notice: this animation has a time period of 0.8 s. In the tutorial case "single_sphere", only 0.01 s is set (tutorial/single_sphere/CFD/system/controlDict line 26: endTime 0.01) just for the validation of successful installation. The user should change the endtime into 0.8 s to get the same result as the above animation.

### Fish settling (corresponds to the tutorial case "multi_sphere_fish"): 
<img src="https://github.com/shenzhihao/CFDEMInterFoamIB/blob/main/animations/fish.gif" width=80% height=80%>  

### Seepage ( not included in tutorial cases):  
<img src="https://github.com/shenzhihao/CFDEMInterFoamIB/blob/main/animations/seepage10mb.gif" width=80% height=80%>  

### Overtopping (not included in tutorial cases):  
<img src="https://github.com/shenzhihao/CFDEMInterFoamIB/blob/main/animations/overtop.gif" width=80% height=80%>  

## Literature about this solver:  
[Shen, Z., Wang, G., Huang, D., & Jin, F. (2022). A resolved CFD-DEM coupling model for modeling two-phase fluids interaction with irregularly shaped particles. Journal of Computational Physics, 448, 110695.](https://www.sciencedirect.com/science/article/pii/S0021999121005908)

## Notice:  
This solver is still on an early release. Questions about installing, compiling, and tutorials cases will be further replied.

## Update logs
2023-05-11:
A case named "multi_sphere_fish" is added in the tutorial. This case can help the users to use the clump DEM model in coupling.
