"""Package M9b: how scipy's two tests behave on ties and all-zero differences."""

import scipy
from scipy.stats import mannwhitneyu, wilcoxon

print("scipy", scipy.__version__)
print("mwu all tied", mannwhitneyu([0, 0, 0], [0, 0, 0], alternative="two-sided"))
try:
    print("wilcoxon all zero", wilcoxon([0, 0, 0], [0, 0, 0]))
except ValueError as e:
    print("wilcoxon all zero raises", e)
print("wilcoxon 3 pairs", wilcoxon([1, 2, 3.5], [0, 0, 0]))
print("wilcoxon 2 pairs", wilcoxon([1, 2], [0, 0]))
