import os

# shortcuts
pj = os.path.join

# configurable stuff
data = config.get('data', 'data')
out = config.get('out', 'out')
trimmomatic_dir = config.get('trimmomatic_dir',
                             'software/trimmomatic-0.36-5/share/trimmomatic')
reads_dir = config.get('reads_dir',
                       './')

# directories
gff = pj(data, 'gff')
roary_dir = pj(out, 'roary')
transcripts_dir = pj(out, 'transcripts')
indexes_dir = pj(out, 'indexes')
counts_dir = pj(out, 'counts')
counts_ref_dir = pj(out, 'counts_reference')
de_dir = pj(out, 'de')
de_ref_dir = pj(out, 'de_reference')

# data files
table = pj(data, 'samples.tsv')
equivalent = pj(data, 'equivalent.tsv')
strains = set([x.rstrip().split('\t')[1]
               for x in open(table)][1:])

# output files
pangenome = pj(roary_dir, 'gene_presence_absence.csv')
transcripts = [pj(transcripts_dir, '%s.fa' % x)
               for x in strains]
indexes = [pj(indexes_dir, '%s.idx' % x)
           for x in strains]
counts = [pj(counts_dir,
             x.rstrip().split('\t')[1],
             x.rstrip().split('\t')[2],
             'abundance.tsv')
          for x in open(table)][1:]
counts_ref = [pj(counts_ref_dir,
                 x.rstrip().split('\t')[1],
                 x.rstrip().split('\t')[2],
                 'abundance.tsv')
              for x in open(table)][1:]
contrasts = sorted({pj(de_dir, '%s.csv' % x)
                    for x in strains
                    if x != 'NT12001'})
contrasts_ref = sorted({pj(de_ref_dir, '%s.csv' % x)
                        for x in strains
                        if x != 'NT12001'})
vst_straight = pj(out, 'vst_straight.tsv')
vst_corrected = pj(out, 'vst_corrected.tsv')

rule pangenome:
  input: gff
  output: pangenome
  params: roary_dir
  threads: 40
  shell:
    'rm -rf {params} && roary -p {threads} -f {params} -s -v -g 100000 {input}/*.gff'

rule transcripts:
  input: transcripts

rule:
  input:
    p=pangenome,
    g=gff,
    e=equivalent
  output: pj(transcripts_dir, '{strain}.fa')
  params: 'NT12001'
  shell:
    'python bin/pangenome2genes.py {input.p} {input.g} {wildcards.strain} --reference {params} --equivalent {input.e} > {output}'

rule indexes:
  input: indexes

rule:
  input: pj(transcripts_dir, '{strain}.fa')
  output: pj(indexes_dir, '{strain}.idx')
  shell: 'kallisto index -i {output} {input}'

rule counts:
  input: counts

rule:
  input:
    index=pj(indexes_dir, '{strain}.idx'),
    rf=table
  output:
    pj(counts_dir, '{strain}', '{replica}', 'abundance.tsv')
  params:
    rdir=reads_dir,
    dir1=pj(counts_dir, '{strain}'),
    dir2=pj(counts_dir, '{strain}', '{replica}'),
    tdir=trimmomatic_dir,
    average=130,
    sd=70
  shell:
    'mkdir -p {params.dir1} && trimmomatic SE {params.rdir}$(awk \'{{if ($2 == "{wildcards.strain}" && $3 == "{wildcards.replica}") print $1}}\' {input.rf}) /dev/stdout ILLUMINACLIP:{params.tdir}/adapters/TruSeq3-SE.fa:2:30:10 | kallisto quant -b 100 -i {input.index} -o {params.dir2} --bias --single -l {params.average} -s {params.sd} /dev/stdin'

rule counts_reference:
  input: counts_ref

rule:
  input:
    index=pj(indexes_dir, 'NT12001.idx'),
    rf=table
  output:
    pj(counts_ref_dir, '{strain}', '{replica}', 'abundance.tsv')
  params:
    rdir=reads_dir,
    dir1=pj(counts_ref_dir, '{strain}'),
    dir2=pj(counts_ref_dir, '{strain}', '{replica}'),
    tdir=trimmomatic_dir,
    average=130,
    sd=70
  shell:
    'mkdir -p {params.dir1} && trimmomatic SE {params.rdir}$(awk \'{{if ($2 == "{wildcards.strain}" && $3 == "{wildcards.replica}") print $1}}\' {input.rf}) /dev/stdout ILLUMINACLIP:{params.tdir}/adapters/TruSeq3-SE.fa:2:30:10 | kallisto quant -b 100 -i {input.index} -o {params.dir2} --bias --single -l {params.average} -s {params.sd} /dev/stdin'

rule de:
  input:
    rf=table,
    c=counts_ref
  output:
    contrasts_ref
  params:
    d=de_ref_dir,
    c=counts_ref_dir
  threads: 40
  shell:
    'Rscript bin/deseq.R {input.rf} {params.c} {params.d} --cores {threads} --pvalue 0.99 --foldchange 0'

rule vst:
  input:
    rf=table,
    c=counts_ref
  output:
    vst_straight,
    vst_corrected
  params:
    counts_ref_dir
  shell:
    'Rscript bin/normalize_counts.R {input.rf} {params} {output}'
