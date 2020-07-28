#ifndef TIDEFAC_GLOBAL_H
#define TIDEFAC_GLOBAL_H

#if defined(_MSC_VER) || defined(WIN64) || defined(_WIN64) ||  \
    defined(__WIN64__) || defined(WIN32) || defined(_WIN32) || \
    defined(__WIN32__) || defined(__NT__)
#define Q_DECL_EXPORT __declspec(dllexport)
#define Q_DECL_IMPORT __declspec(dllimport)
#elif defined(SWIG)
#define Q_DECL_EXPORT
#define Q_DECL_IMPORT
#else
#define Q_DECL_EXPORT __attribute__((visibility("default")))
#define Q_DECL_IMPORT __attribute__((visibility("default")))
#endif

#if defined(TIDEFAC_LIBRARY)
#define TIDEFAC_EXPORT Q_DECL_EXPORT
#else
#define TIDEFAC_EXPORT Q_DECL_IMPORT
#endif

#endif  // TIDEFAC_GLOBAL_H
