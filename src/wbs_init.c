#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .C calls */
extern void bs_rec_wrapper(void *, void *, void *);
extern void wbs_int_rec_wrapper(void *, void *, void *, void *, void *);
extern void wbs_rec_wrapper(void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"bs_rec_wrapper",      (DL_FUNC) &bs_rec_wrapper,      3},
    {"wbs_int_rec_wrapper", (DL_FUNC) &wbs_int_rec_wrapper, 5},
    {"wbs_rec_wrapper",     (DL_FUNC) &wbs_rec_wrapper,     5},
    {NULL, NULL, 0}
};

void R_init_wbs(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
