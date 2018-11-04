/*
 * This file contains inline utility functions and function declarations needed
 * for SPH
 * 
 * Based on the physics modules implemented by James Wadsley in GASOLINE.
 * Implemented by Tom Quinn in ChaNGa and refactored by Isaac Backus
 * Oct 14, 2016
 */
#ifndef SPHUTILS_H
#define SPHUTILS_H

// Include Sph.h for possible macro definitions in there
#include "Sph.h"

/**
 * @brief presPdv can be used as a simple means to switch Pressure weighting
 * schemes for PdV calculations.
 * @return a
 */
inline double presPdv(const double a, const double b){
    return a;
}
/**
 * @brief presAcc can be used as a simple means to switch Pressure weighting for
 * schemes for acceleration calculation.
 * @return a + b
 */
inline double presAcc(const double a, const double b){
    return a + b;
}

/**
 * @brief switchCombine is used to allow weighting of the BalsaraSwitch of 
 * two particles.
 */
inline double switchCombine(GravityParticle *a, GravityParticle *b) {
#ifdef CULLENALPHA
    return 1.0;
#else
    return 0.5*(a->BalsaraSwitch()+b->BalsaraSwitch());
#endif
}

/**
 * @brief rhoDivv can be used for switching between using rtforce and not
 * @return b if rtforce is compiled, otherwise return a
 */
inline double rhoDivv(const double a, const double b) {
#ifdef RTFORCE
    return b;
#else
    return a;
#endif
}

/**
 * @brief varAlpha allows for variable SPH viscosity parameter alpha. Behavior
 * is controlled by compile-time flag VARALPHA.
 */
inline double varAlpha(const double alpha, const GravityParticle *a, 
                       const GravityParticle *b) {
#if defined(VARALPHA) || defined(CULLENALPHA)
    return alpha*0.5*(a->CullenAlpha() + b->CullenAlpha());
#else
    return alpha;
#endif //VARALPHA
}

/**
 * @brief varBeta allows for variable SPH viscosity beta.  Behavior is controlled
 * by the compile-time flag VARALPHA
 */
inline double varBeta(const double beta, GravityParticle *a, 
                      const GravityParticle *b) {
#if defined(VARALPHA) || defined(CULLENALPHA)
    return beta*0.5*(a->CullenAlpha() + b->CullenAlpha());
#else
    return beta;
#endif //VARALPHA
}

/**
 * @brief diffusionLimitTest performs a test which can be used to prevent the
 * diffusion of supernova feedback energy.  Its behavior is controlled by
 * FEEDBACKDIFFLIMIT
 * @param diffSum The possibly weighted sum of interacting particles diff()
 * @param dTime Current time
 * @param a Particle 1
 * @param b Neighbor particle
 */
inline bool diffusionLimitTest(const double diffSum, const double dTime, 
                               GravityParticle *a, GravityParticle *b) {
#ifdef FEEDBACKDIFFLIMIT
    return (diffSum == 0 \
            || dTime < a->fTimeCoolIsOffUntil() \
            || dTime < b->fTimeCoolIsOffUntil());
#else
    return (diffSum == 0);
#endif
}

inline double massDiffFac(const GravityParticle *p) {
    return 1.0;
}

typedef struct PressSmoothUpdateStruct {
    double dvdotdr;
    double visc;
    double aFac;
    Vector3D<double> dx;
    #ifdef DIFFUSION
            #ifndef NODIFFUSIONTHERMAL /* compile-time flag */
                double diffu;
            #endif
        double diffMetals;
        double diffMetalsOxygen;
        double diffMetalsIron;
    #endif
} PressSmoothUpdate;

typedef struct PressSmoothParticleStruct {
    double rNorm;
    double PoverRho2;
    double PoverRho2f;
} PressSmoothParticle;

void updateParticle(GravityParticle *a, GravityParticle *b, 
                    PressSmoothUpdate *params, PressSmoothParticle *aParams, 
                    PressSmoothParticle *bParams, int sign);

#endif // SPHUTILS_H
