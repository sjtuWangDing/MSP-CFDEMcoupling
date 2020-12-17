#ifndef __OF_VERSION_H__
#define __OF_VERSION_H__

#if (MSP_CFDEM_WM_PROJECT_VERSION == 50)
#define version50
#elif (MSP_CFDEM_WM_PROJECT_VERSION == 40)
#define version40
#elif (MSP_CFDEM_WM_PROJECT_VERSION == 24)
#define version24
#elif (MSP_CFDEM_WM_PROJECT_VERSION == 30)
#define version30
#elif (MSP_CFDEM_WM_PROJECT_VERSION == 132)
#define versionExt32
#elif (MSP_CFDEM_WM_PROJECT_VERSION == 1606)
#define versionv1606plus
#elif (MSP_CFDEM_WM_PROJECT_VERSION == 1612)
#define versionv1612plus
#elif (MSP_CFDEM_WM_PROJECT_VERSION == 1706)
#define versionv1706
#endif

// define anisotropicRotation cloud models
// #define anisotropicRotation

// features of 16ext work also in extend 3.2
#if defined(versionExt32)
#define version16ext
#endif

// features of 4.0 work also in 5.0
#if defined(version50)
#define version40
#endif

// features of 3.0 work also in 4.0
#if defined(version40)
#define version30
#endif

// features of v1706 work also in v1612+
#if defined(versionv1706)
#define versionv1612plus
#endif

#if defined(versionv1612plus)
#define version40
#endif

// features of v1606+ work also in v1612+
#if defined(versionv1612plus)
#define versionv1606plus
#endif

// features of 3.0 work also in v1606+
#if defined(versionv1606plus)
#define version30
#endif

// features of 2.4Dev work also in 3.0
#if defined(version30)
#define version24Dev
#endif

// basically use 24x settings + some dev features (e.g. new turbulence model structure)
#if defined(version24Dev)
#define version24
#endif

// features of 2.4 work also in 2.3
#if defined(version24)
#define version23
#endif

// features of 2.1 work also in 2.3
#if defined(version23)
#define version21
#define version221
#endif

// features of 2.1 work also in 2.2
#if defined(version22)
#define version21
#define version221
#endif

#endif  // __OF_VERSION_H__
