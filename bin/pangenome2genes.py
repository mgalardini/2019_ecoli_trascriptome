#!/usr/bin/env python


import os
import sys
import argparse
import pandas as pd
from BCBio import GFF


def get_options():
    # create the top-level parser
    description = "Extract putative transcripts from a series of gff files"
    parser = argparse.ArgumentParser(description=description)
    
    parser.add_argument('pangenome',
                        help='Roary\'s gene presence absence file')
    parser.add_argument('gffdir',
                        help='GFF files directory (format: COLUMN_ID.gff)')
    parser.add_argument('priority',
                        help='Priority strain (i.e. Roary\'s column)')
    
    parser.add_argument('--reference',
                        default=None,
                        help='Only consider genes from a specific organism [Default: full pangenome]')
    parser.add_argument('--equivalent',
                        default=None,
                        help='Equivalences in roary\'s columns, in tabular format '
                             '(column -> equivalent_strain). '
                             'Useful if two or more GFF files are identical '
                             'when running roary')

    return parser.parse_args()


if __name__ == "__main__":    
    options = get_options()

    # Load roary
    roary = pd.read_table(options.pangenome,
                          sep=',',
                          low_memory=False)
    # Set index (group name)
    roary.set_index('Gene', inplace=True)
    # Drop the other info columns
    roary.drop(list(roary.columns[:13]), axis=1, inplace=True) 

    # Focus on a "reference" strain only?
    if options.reference is not None:
        roary = roary.loc[roary[options.reference].dropna().index]
        roary.index = roary[options.reference]
        roary.index.name = 'Gene'
    
    # Consider equivalent strains
    equi = {}
    for l in open(options.equivalent):
        orig, other = l.rstrip().split('\t')
        if other not in roary.columns:
            roary[other] = roary[orig].copy()

    # Prepare all the GFF files parsers
    parsers = {}
    for column in roary.columns:
        parsers[column] = GFF.parse(os.path.join(options.gffdir,
                                                 '%s.gff' % column))
    # Parse the "priority" strain
    genes = {}
    gnames = roary[
            options.priority
            ].dropna().reset_index().set_index(
                    options.priority
                    )['Gene'].to_dict()
    # There might still be some in-paralogs there
    add = {}
    for k,v in gnames.items():
        if '\t' in k:
            for k1 in k.split('\t'):
                add[k1] = v
        if '___' in k:
            for k1 in k.split('\t'):
                add[k1.split('___')[0]] = v
    for k, v in add.items():
        gnames[k] = v
    genes[options.priority] = {}
    unnamed = 1
    for s in parsers[options.priority]:
        for f in filter(lambda x: x.type == 'gene',
                        s.features):
            # Name
            # Try first with the locus_tag
            try:
                name = f.qualifiers['locus_tag'][0]
            except:
                name = 'unnamed_%d' % unnamed
                unnamed += 1

            name = gnames.get(name, name)

            # Sequence
            seq = str(f.extract(s).seq)

            genes[options.priority][name] = seq

    # Go gene by gene
    for gene in roary.index:
        agenes = roary.loc[gene].dropna()
        if options.priority in agenes.index:
            print('>%s\n%s' % (gene, genes[options.priority][gene]))
        else:
            # Pick the first column, parse it
            strain = list(agenes.index)[0]
            gnames = roary[
                    strain
                    ].dropna().reset_index().set_index(
                            strain
                            )['Gene'].to_dict()
            # There might still be some in-paralogs there
            add = {}
            for k,v in gnames.items():
                if '\t' in k:
                    for k1 in k.split('\t'):
                        add[k1] = v
                if '___' in k:
                    for k1 in k.split('\t'):
                        add[k1.split('___')[0]] = v
            for k, v in add.items():
                gnames[k] = v
            if strain not in genes or len(genes[strain]) == 0:
                # parse
                genes[strain] = {}
                for s in parsers[strain]:
                    for f in filter(lambda x: x.type == 'gene',
                                    s.features):
                        # Name
                        # Try first with the locus_tag
                        try:
                            name = f.qualifiers['locus_tag'][0]
                        except:
                            name = 'unnamed_%d' % unnamed
                            unnamed += 1

                        name = gnames.get(name, name)

                        # Sequence
                        seq = str(f.extract(s).seq)

                        genes[strain][name] = seq
            print('>%s\n%s' % (gene, genes[strain][gene]))
