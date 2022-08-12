In this version, only case 1 is completed as final version. As
previously, kd1, kd2 and kex parameters were appropriately fixed. 
The number of outer iterations are increased to 200, and the 
number of inner iterations are decreased to 10. 

In Test Run 1, parameters are randomized within 0.9 to 1.1 times 
of P0. The parameters are bounded to be within P/5 and P*5. All 
parameters and model fits are saved under Saved_Parameters1. 

In Test Run 2 (final), parameters are randomized within 0.8 
to 1.2 times of P0. The parameters are bounded to be within 
P/10 and P*10. All parameters and model fits are saved under 
Saved_Parameters2. 

The idea was to reduce the spread of the parameter space 
(i.e. to reduce SD) over the outer iterations.

Also AIC calculations are introduced and the code is little
restructured to save calculated parameters in the diary with
nice format so that things can be easily integrated into
the manuscript without ambiguity.

The initial guesses for the regulatory parameters are such 
that there is no effect initially with the regulatory terms.
Initial guesses are as follows. 
M1: p0 = [6E2,9E2,1E-1,6E2,9E3,1E-1,1E-6,3E-5,1.5E0];
M2: p0 = [6E2,9E2,1E-1,1E5,6E2,9E3,1E-1,1E5,1E-6,3E-5,1.5E0];
M3: p0 = [6E2,9E2,1E-1,1E1,6E2,9E3,1E-1,1E1,1E-6,3E-5,1.5E0];
M4: p0 = [6E2,9E2,1E-1,5E4,6E2,9E3,1E-1,2.5E4,1E-6,3E-5,1.5E0];
M5: p0 = [6E2,9E2,1E-1,5E4,1E5,6E2,9E3,1E-1,2.5E4,1E5,1E-6,3E-5,1.5E0];
M6: p0 = [6E2,9E2,1E-1,5E4,1E1,6E2,9E3,1E-1,2.5E4,1E1,1E-6,3E-5,1.5E0];


This is the final version for PNAS paper revision.
