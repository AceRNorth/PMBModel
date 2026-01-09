
#ifndef GENERALMETAPOP_CONSTANTS_H
#define GENERALMETAPOP_CONSTANTS_H

/**
 * Namespace for constants used throughout the program. 
 */
namespace constants 
{
    const int max_dev = 20; /**< Juvenile development time (egg to adult) expressed as days left till eclosion (eclosion on day 0). */ 
    const int num_gen = 3; /**< The number of different genotypes in the mosquito population. For 6 genotypes these represent WW, WD, DD, WR, RR and DR (in this order), where W is the wild type allele, D is the drive type allele and R is the non-functional resistance allele. */ 
    const double pi = 3.14159265; /**< Numerical constant PI (to 8 d.p.). */
}

/** 
 * An enum for the allowed boundary type values in the model. 
 */
enum BoundaryType{
    Toroid, /**< Toroidal (or periodic) boundary conditions. @see ToroidalBoundaryStrategy */
    Edge /**< Edge boundary conditions. @see EdgeBoundaryStrategy */
}; 

/** 
 * An enum for the allowed dispersal type values in the model. 
 */
enum DispersalType{
    DistanceKernel, /**< Distance kernel dispersal type. @see DistanceKernelDispersal */
    Radial /**< Radial dispersal type. @see RadialDispersal */
};

#endif //GENERALMETAPOP_CONSTANTS_H
