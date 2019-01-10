#!/usr/bin/env python


import os
import shutil
import argparse


def get_options():
    description = 'Copy input GFF files to data directory'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('samples',
                        help='Samples file')
    parser.add_argument('gff',
                        help='Gff directory')
    parser.add_argument('out',
                        help='Output directory')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    strains = set()
    b = True
    for l in open(options.samples):
        if b:
            b = False
            continue
        strain = l.rstrip().split('\t')[1]
        strains.add(strain)

    for strain in sorted(strains):
        if strain == 'NT12001':
            gff = 'genome.gff'
        else:
            files = [x for x in os.listdir(options.gff)
                     if x.startswith(strain)]
            if len(files) == 0:
                print('Genome %s not found!' % strain)
                continue
            gff = sorted(files)[0]
        shutil.copy(os.path.join(options.gff, gff),
                    os.path.join(options.out, '%s.gff' % strain))
