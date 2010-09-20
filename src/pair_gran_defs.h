    //use #define here, so we dont allocate new properties unnecessarily
    //Atom types start from 1, not zero (this is ensured by fix/pour)

#ifdef LMP_GRAN_DEFS_DEFINE
    #define itype (atom->type[ip])
    #define jtype (atom->type[jp])
    #define ri (atom->radius[ip])
    #define rj (atom->radius[jp])
    #define reff_wall (atom->radius[ip])
    #define sqrtFiveOverSix 0.91287092917527685576161630466800355658790782499663875
#elif defined LMP_GRAN_DEFS_UNDEFINE
    #undef ri
    #undef rj
    #undef reff_wall
    #undef itype
    #undef jtype
#endif
