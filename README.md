# Lorenz-System

*lorenz.m* contains the instructions and algorithm to find and compute the aproximate inestable periodic orbit of the Lorenz system for r = 18 and r = 21 (To compute different r's, define r and a good starting point x0).

Complementary functions:

  1. *P.m*
  
   Given the vectors *X* and *T* that contain the corresponding points and times obtained with *ode45* of an orbit, a starting point *x* and a plane *&Sigma;:{z = h}*, it finds the exact point in *&Sigma;* that corresponds to compute a return from the starting point *x* and its orbital period.
   
   2. *Q.m*

   Given a starting point *x* and a plane *&Sigma;:{z = h}*, it finds the euclidean distance between the starting point and the one in *&Sigma;* obtained in *P.m*.


The main script *lorenz.m* is divided into different parts:

  1. 
