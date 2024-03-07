# AEP

This repository contains the **original C++ package** to use the AEP algorithm described in the paper

**Arbenz, P., Embrechts, P., and G. Puccetti (2011). The AEP algorithm for the fast computation of the distribution of the sum of dependent random variables  available at https://doi.org/10.3150/10-BEJ284**

The AEP algorithm numerically calculates the distribution function and the Value-at-Risk of the sum of a number of d positive random variables, when the joint distribution H of the vector (X1,...,Xd) is given. 

**AEP.C is the C++ code for AEP**. Note that its not particularly beautiful code, 
and it will probably be difficult to understand the ideas behind all the pointers and 
vectors. however, the good thing is that you only have to adapt the dimension (DIM) and 
the JDF (specify the marginals and copula). note that the copula is required to be 
analytic.

In the folder _executable_, one can find the **executable files** for Linux, Mac and Windows michines.
AEP-Df* computes the distribution function, AEP-VaR* the Value-at-Risk

**PLEASE READ THE USER GUIDE AVAILABLE IN PDF FORMAT**

Notice that a new **Python** library including an implementation of the AEP algorithm 
has been provided in https://github.com/gpitt71/gemact-code
