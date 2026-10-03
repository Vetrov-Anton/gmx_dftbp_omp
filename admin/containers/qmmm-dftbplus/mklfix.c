/*
 * Override of the CPU vendor check of Intel MKL, decided at run time.
 *
 * MKL takes its optimized code paths (AVX2, AVX-512) only if its internal function
 * mkl_serv_intel_cpu_true() returns 1, which it does for Intel CPUs only; on AMD CPUs MKL
 * falls back to slower generic paths (MKL_DEBUG_CPU_TYPE has been ignored since MKL 2020.1).
 * This library defines the same symbol and is loaded ahead of MKL (it is a DT_NEEDED of gmx
 * and dftb+, see gmx-dftbplus-mkl.def), so MKL calls this function instead:
 *   - on AMD (and Hygon) CPUs, found with CPUID, it returns 1;
 *   - on any other CPU it returns what the function of MKL returns (looked up with RTLD_NEXT),
 *     i.e. MKL behaves as without this library;
 *   - MKL_VENDOR_OVERRIDE=off leaves MKL alone on any CPU, MKL_VENDOR_OVERRIDE=on forces 1.
 * Built with -DMKLFIX_ALWAYS (CPU_VENDOR=amd), it returns 1 unconditionally.
 */
#define _GNU_SOURCE
#include <cpuid.h>
#include <dlfcn.h>
#include <stdlib.h>
#include <string.h>

static int cpuIsAmd(void)
{
    unsigned int eax, ebx, ecx, edx;
    char         vendor[13];
    if (!__get_cpuid(0, &eax, &ebx, &ecx, &edx))
    {
        return 0;
    }
    memcpy(vendor, &ebx, 4);
    memcpy(vendor + 4, &edx, 4);
    memcpy(vendor + 8, &ecx, 4);
    vendor[12] = '\0';
    return strcmp(vendor, "AuthenticAMD") == 0 || strcmp(vendor, "HygonGenuine") == 0;
}

static int mklOwnAnswer(void)
{
    int (*mklFunction)(void) = (int (*)(void))dlsym(RTLD_NEXT, "mkl_serv_intel_cpu_true");
    return mklFunction != NULL ? mklFunction() : 0;
}

int mkl_serv_intel_cpu_true(void)
{
    static int answer = -1; /* decided once; a race only computes the same value twice */
    if (answer < 0)
    {
        const char* env = getenv("MKL_VENDOR_OVERRIDE");
#ifdef MKLFIX_ALWAYS
        int force = 1;
#else
        int force = cpuIsAmd();
#endif
        if (env != NULL && strcmp(env, "off") == 0)
        {
            force = 0;
        }
        else if (env != NULL && strcmp(env, "on") == 0)
        {
            force = 1;
        }
        answer = force ? 1 : mklOwnAnswer();
    }
    return answer;
}
