# Potentials

A plot of the Hamiltonian and potential energy surfaces overlaid withe the Hamiltonian vector field (**Figure 2**) was generated using:

```
PotentialField2(2,'steps',20,'startPoint',[1,1],'npoints',1000,'step_size',1e-3,'save',true)
```

A plot of a series of tranjectories tilted out of the plane, together with a top-down projection, (**Figure 3**) was generated using:

```
PotentialField8(2,'steps',5,'npoints',1000,'step_size',1e-3,'save',true)
```

A plot of a series of three-dimensional trajectories (**Figure 4**A) was generated using

```
PotentialField7(3,'steps',5,'npoints',3000,'step_size',1e-3,'save',true)
```

and a plot of the a series of projected three-dimensional trajectories (**Figure 4**B) was generated using

```
PotentialField8(3,'steps',5,'npoints',3000,'step_size',1e-3,'save',true)
```

Plots of the distortion induced by a third reference state in the plane (**Figure 5**) were generated for trajectories started at different points were generated using (A)

```
PotentialField3(2,'steps',20,'startPoint',[0.5,0.5],'extraStates',[1,1],'npoints',1000,'step_size',1e-3,'magnification',100,'save',true)
```

and (B)

```
PotentialField3(2,'steps',20,'startPoint',[1,1],'extraStates',[1,1],'npoints',1000,'step_size',1e-3,'magnification',100,'save',true)
```

A plot of the potential energy surface (**Figure S6**A) and conjugate potential surface (**Figure S6**B) overlaid with the potential field using

```
PotentialField10(2,'steps',20,'startPoint',[1/(2)+0.1,1/(2)],'npoints',1200,'step_size',1e-3,'save',true)
```
