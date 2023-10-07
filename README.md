# Run
```
root math_problem.cxx
```

# Description
 
This is the implementation of simple physics problem:
 * mass 1, originally at rest, v1_0=0
 * mass 2, originally at speed -v2_0 (from right to left)
Placed as:
 
```
// ||| 
// |||____|m1|____<--|m2|___
// 
// ||| = a wall that reflects m1 without energy loss
```

# Statement of the problem

The number of collisions needed to reflect back m2 tends asymtotically to pi/2*sqrt(m2/m1)


```
  m2		  2*collisions/sqrt(m2)		 Real Time (s)
1.00e+00	6.0000000			0.0000019
1.00e+02	3.4000000			0.0000021
1.00e+04	3.1800000			0.0000100
1.00e+06	3.1440000			0.0000861
1.00e+08	3.1418000			0.0008452
1.00e+12	3.1415960			0.0855410
1.00e+14	3.1415930			0.8539388
```

# Analytical solution

I found it later in [youtube](https://youtu.be/jsYwFizhncE)
TLDR: solve the problem in the phase space
