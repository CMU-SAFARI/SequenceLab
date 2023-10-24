import sys
import numpy as np
BOTTOM_SPLIT = 90


if __name__ == "__main__":
    if len(sys.argv) >= 3:
        filename = sys.argv[1]
        tss_filename = sys.argv[2]
    else:
        print("no filename for lengths and/or tss provided provided")
        exit(-1)

    lengths = np.loadtxt(filename, dtype=int)

    cutoff = np.percentile(lengths, BOTTOM_SPLIT)

    top = np.argwhere(lengths > cutoff)
    print(cutoff)

    bottom_filename = tss_filename.rsplit(".", maxsplit=1)[0] + "_bottom.tss"
    top_filename = tss_filename.rsplit(".", maxsplit=1)[0] + "_top.tss"

    with open(tss_filename, "r") as tss:
        with open(bottom_filename, "w") as bottom_file:
            with open(top_filename, "w") as top_file:
                idx = 0
                for line in tss:
                    if idx in top:
                        top_file.write(line)
                    else:
                        bottom_file.write(line)
                    idx += 1
    print("Done!")
