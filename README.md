# Potentials

A plot of the Hamiltonian and potential energy surfaces overlaid withe the Hamiltonian vector field (**Figure 2**) was generated using:

```
PotentialField2(2,'steps',20,'startPoint',[1,1],'npoints',1000,'step_size',1e-3,'save',true)
```

A plot of a series of tranjectories tilted out of the plane, together with a top-down projection, (**Figure 3**) was generated using:

```
PotentialField8(2,'steps',5,'npoints',1000,'step_size',1e-3,'save',true)
```

and a plot of the series of three-dimensional suraces for increasing values of the inter-state separation distance c (**Figure 4**) was generated using

```
PotentialField11(3,1,'steps',100,'npoints',3000,'step_size',1e-3,'save',true)
PotentialField11(3,0.5,'steps',100,'npoints',3000,'step_size',1e-3,'save',true)
PotentialField11(3,0.25,'steps',100,'npoints',3000,'step_size',1e-3,'save',true)
PotentialField11(3,0.2,'steps',100,'npoints',3000,'step_size',1e-3,'save',true)
PotentialField11(3,0.125,'steps',100,'npoints',3000,'step_size',1e-3,'save',true)
PotentialField11(3,0.063,'steps',100,'npoints',3000,'step_size',1e-3,'save',true)
```
Calculation of the inverse of the combination matrix B was perfomed using

```
combination_matrix(3)
combination_matrix(4)
```
for use in the Wolfram Notebook <kbd>Moment of inertia matrix.nb</kbd>, wherein the n × n matrices for the three- and four-state systems are computed symbollically.
