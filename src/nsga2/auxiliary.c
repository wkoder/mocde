/* Some utility functions (not part of the algorithm) */

# include "global.h"

#ifdef NSGA2_IMPL

# include <stdio.h>
# include <stdlib.h>
# include <math.h>


/* Function to return the maximum of two variables */
double maximum (double a, double b)
{
    if (a>b)
    {
        return(a);
    }
    return (b);
}

/* Function to return the minimum of two variables */
double minimum (double a, double b)
{
    if (a<b)
    {
        return (a);
    }
    return (b);
}

#endif
