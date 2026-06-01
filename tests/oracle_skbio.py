"""scikit-bio ANOSIM oracle. Reads a square distance-matrix TSV and a grouping
TSV, prints `field<TAB>value` lines matching rsomics-anosim's output. Used by
tests/compat.rs for the live differential."""
import sys
import numpy as np
from skbio import DistanceMatrix
from skbio.stats.distance import anosim


def read_dm(path):
    with open(path) as f:
        header = None
        ids = []
        rows = []
        for line in f:
            if line.strip() == "" or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                ids = [p.strip() for p in parts[1:]]
                continue
            rows.append([float(x) for x in parts[1:]])
    return DistanceMatrix(np.array(rows), ids)


def read_grouping(path, ids):
    m = {}
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s == "" or s.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            m[parts[0].strip()] = parts[1].strip()
    return [m[i] for i in ids]


def main():
    dm_path, grp_path = sys.argv[1], sys.argv[2]
    perms = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 42
    dm = read_dm(dm_path)
    grouping = read_grouping(grp_path, dm.ids)
    res = anosim(dm, grouping, permutations=perms, seed=seed)
    print("method name\t%s" % res["method name"])
    print("test statistic name\t%s" % res["test statistic name"])
    print("sample size\t%d" % res["sample size"])
    print("number of groups\t%d" % res["number of groups"])
    print("test statistic\t%r" % float(res["test statistic"]))
    print("p-value\t%r" % float(res["p-value"]))
    print("number of permutations\t%d" % res["number of permutations"])


if __name__ == "__main__":
    main()
