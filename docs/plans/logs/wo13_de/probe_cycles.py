"""Package wo13-de: do pydeseq2's objects outlive their last name until the cycle collector runs?

If a DeseqDataSet sits in a reference cycle, its layers stay allocated after ``fit`` returns
until ``gc`` happens to run: a large numpy array does not count towards gc's thresholds.

Run: bash env.sh python docs/plans/logs/wo13_de/probe_cycles.py
"""

import gc
import weakref

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

gc.disable()
rng = np.random.default_rng(0)
samples = [f"s{i}" for i in range(12)]
genes = [f"g{i}" for i in range(50)]
counts = pd.DataFrame(rng.poisson(20, (12, 50)), index=samples, columns=genes)
metadata = pd.DataFrame({"group": ["a"] * 6 + ["b"] * 6}, index=samples)
dds = DeseqDataSet(counts=counts, metadata=metadata, design="~group", quiet=True)
dds.deseq2()
stats = DeseqStats(dds, contrast=["group", "b", "a"], quiet=True)
stats.summary()
alive = weakref.ref(dds), weakref.ref(stats)
del dds, stats
print("after del: dds alive", alive[0]() is not None, "; stats alive", alive[1]() is not None)
gc.collect()
print("after gc.collect: dds alive", alive[0]() is not None, "; stats", alive[1]() is not None)
