# Lorenz-System

*lorenz.m* contains the instructions and algorithm to find and compute the aproximate inestable periodic orbit of the Lorenz system for r = 18 and r = 21 (To compute a different r, define r and a good starting point x0).

Complementary functions:

  1. *P.m*
  
   Given the vectors *X* and *T* that contain the corresponding points and times obtained with *ode45* of an orbit, a starting point *x* and a plane *&Sigma;:{z = h}*, it finds the exact point in *&Sigma;* that corresponds to compute a return from the starting point *x* and its orbital period.
   
   2. *Q.m*

   Given a starting point *x* and a plane *&Sigma;:{z = h}*, it finds the euclidean distance between the starting point and the one in *&Sigma;* obtained in *P.m*.


The main script *lorenz.m* is divided into different parts:

  1. It finds for which *r* the eigenvalues of *Df* corresponding to the Lorenz system go from being real to complex.
  2. It computes the pitchfork bifurcation in *r = 1*.
  3. Finds the inestable periodic orbit for *r = 18*.
  4. Finds the inestable periodic orbit for *r = 21*.


Figures legend:

1. Orbit values of *x*, *y* and *z* along time for x0 defined above for *r = 18*.
2. 3DPlot of the orbit starting at x0 for *r=18*. The plane *&Sigma;* is computed in blue (*r = 18*).
3. XY representation of the inestable periodic orbit (*r = 18* and *r = 21*)
4. Orbit values of *x*, *y* and *z* along time for a point in the inestable periodic orbit (*r = 18*).
5. 3DPlot of the inestable periodic orbit (green) and the orbit starting at x1 defined above (*r = 18*).
6. Orbit values of *x*, *y* and *z* along time for x0 defined above for *r = 21*.
7. 3DPlot of the orbit starting at x0 for *r=18*. The plane *&Sigma;* is computed in blue (*r = 21*).
8. Orbit values of *x*, *y* and *z* along time for a point in the inestable periodic orbit (*r = 21*).
9. 3DPlot of the inestable periodic orbit (green) and the orbit starting at x1 defined above (*r = 21*).
10. Points used for polynomial interpolation in P function.
11. Representation of the Pitchfork bifurcacion.
12. Plot of the Hopf subcritical bifurcation (*r aprox 24.73...*) -> x coordinate
13. Plot of the Hopf subcritical bifurcation (*r aprox 24.73...*) -> y coordinate
