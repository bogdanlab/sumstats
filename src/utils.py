import scipy.stats as st
import numpy as np

def ptoz(pval, sign):
    z = np.fabs(st.norm.ppf(pval/2.0))
    if sign == '+':
        return z
    return -1.0*z
