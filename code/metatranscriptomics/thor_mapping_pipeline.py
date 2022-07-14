#!/bin/env python

import glob
import os
import sys

threads = 74

ref = ['fna/Bc.fna', 'fna/Fj.fna', 'fna/Pk.fna']

def check_dir(d):
    if not os.path.isdir(d):
        os.makedirs(d)

def compress(f):
    pigz = ' '.join(['pigz',
                     '-p', str(threads),
                     '-q',
                     f])
    os.system(pigz)

for fq1 in glob.glob("fq/*_R1_001.fastq.gz"):
    fq2 = fq1.replace('_R1_001.fastq.gz', '_R2_001.fastq.gz')
    pref = fq1.split('/')[-1].split('_')[0]
    print('Processing '+pref+':')
    
    ## fastp
    print("\tfastp")
    check_dir('fastp')
    fp1, fp2 = 'fastp/'+pref+'.fp1', 'fastp/'+pref+'.fp2'
    if not os.path.isfile(fp1):
        t = threads
        if t > 16:
            t = 16
        fpcmd = ' '.join(['fastp',
                          '--in1', fq1,
                          '--out1', fp1,
                          '--in2', fq2,
                          '--out2', fp2,
                          '-h', 'fastp/'+pref+'.html',
                          '-j', 'fastp/'+pref+'.json',
                          '-w', str(t)])
        os.system(fpcmd)

    ## bbsplit
    print("\tbbsplit")
    check_dir('bbsplit')
    if not os.path.isfile('bbsplit/'+pref+'_Bc.fq.gz'):
        bbcmd = ' '.join(['bbsplit.sh',
                          'in1='+fp1,
                          'in2='+fp2,
                          'ref='+','.join(ref),
                          'basename=bbsplit/'+pref+'_%.fq',
                          't='+str(threads),
                          'int=t', 'ambig2=toss', 'vslow'])
        os.system(bbcmd)
        for q in glob.glob('bbsplit/'+pref+"*.fq"):
            compress(q)

    ## bowtie indicies
    print("\tbowtie2 indexing")
    check_dir('bt2i')
    for r in ref:
        rp = r.split('/')[-1][:-4]
        if not os.path.isfile('bt2i/'+rp+'.1.bt2'):
            bti = ' '.join(['bowtie2-build',
                            '-f', r,
                            'bt2i/'+rp])
            os.system(bti)
    
    ## mapping and counting
    check_dir('bt2')
    check_dir('htseq')
    for gz in glob.glob("bbsplit/"+pref+"*.fq.gz"):
        ## bt2
        print("\tbowtie2 mapping")
        r = gz.split('/')[-1][-8:-6]
        sam = 'bt2/'+'_'.join([pref, r])+'.sam'
        bam = sam.replace('.sam', '.bam')
        if not os.path.isfile(bam):
            index = 'bt2i/'+r
            bt2 = ' '.join(['bowtie2',
                            '-p', str(threads),
                            '--sensitive',
                            '-x', index,
                            '--interleaved', gz,
                            '-S', sam])
            os.system(bt2)
            os.system('samtools sort '+sam+' > '+bam)
            os.system('rm '+sam)
            os.system('samtools index '+bam)
        ## htseq
        print("\thtseq counting")
        htout = 'htseq/'+'_'.join([pref, r])+'.txt'
        if not os.path.isfile(htout):
            gff = 'as_parse/'+r+'.gff'
            htcmd = ' '.join(['htseq-count',
                              '--stranded', 'reverse',
                              '-m', 'intersection-strict',
                              '-t', 'CDS',
                              '-i', 'ID',
                              '-n', str(threads),
                              bam,
                              gff,
                              '>', htout])
            os.system(htcmd)
