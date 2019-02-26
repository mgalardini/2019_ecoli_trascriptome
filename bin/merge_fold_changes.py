#!/usr/bin/env python


import os
import sys
import pandas as pd


def get_options():
    import argparse

    description = 'Merge fold changes in a single file'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('input',
                        nargs='+',
                        help='DESeq2 output table')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    m = None
    for f in options.input:
        df = pd.read_csv(f,
                         index_col=0)
        sname = os.path.split(f)[-1].split('.')[0]
        df = df.rename(columns={'log2FoldChange': sname})
        if m is None:
            m = df[sname].to_frame()
        else:
            m = m.join(df[sname].to_frame())

    m.to_csv(sys.stdout,
             sep='\t')
