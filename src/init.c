#include "header.h"
static const R_CallMethodDef CallEntries[] = {
    {"C_get_log_perms", (DL_FUNC) &C_get_log_perms, 7},
    {"C_get_log_perms_bioassay", (DL_FUNC) &C_get_log_perms_bioassay, 8},
    {"C_get_log_ML", (DL_FUNC) &C_get_log_ML, 4},
    {"C_get_log_ML_bioassay", (DL_FUNC) &C_get_log_ML_bioassay, 7},
    {"C_wrapper_log_sum_exp", (DL_FUNC) &C_wrapper_log_sum_exp, 2},
    {NULL, NULL, 0}
};

void R_init_perms(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}