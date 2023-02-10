import numpy as np
import gzip
import itertools


def parseVCF(filename):
    '''read genotype info from VCF format
    Note: this version could only read VCF file all at once!
    '''
    genotype_arr = []
    # with gzip.open(filename, 'rt') as input: this should be adapted in read situation
    with open(filename, 'rt') as input:
        for line in input.readlines():
            if '#' in line:
                continue
            # print( line.strip().split('\t')[9:] )
            # genotype_arr.append([i.split(':')[0] for i in line.strip().split('\t')[9:]])
            genotype_arr.extend([i.split(':')[0].split('/') for i in line.strip().split('\t')[9:]])
    genotype_arr_flatten = list(itertools.chain(*genotype_arr))
    return genotype_arr_flatten

genotype_arr = parseVCF('../testdataset/827657.CEU.biallelicSNPs.chr1.vcf.gz')
# print(genotype_arr)
# print(len(genotype_arr))

num_sequences = len(genotype_arr)
num_sites = 1  # this variable could be passed from command line
segregating_site=1 # this value represent how many SNP in the targeted region (this region will have non-SNP)

### calculation of seperate variables
a1 = (1 / np.arange(1, num_sequences)).sum()
a2 = (1 / np.arange(1, num_sequences)**2).sum()
b1 = (num_sequences+1) / (3*(num_sequences-1))
b2 = 2*(num_sequences**2 + num_sequences + 3)/ (9*num_sequences*(num_sequences-1))
c1 = b1 - (1/a1)
c2 = b2 - (num_sequences+2)/(a1*num_sequences) + a2/(a1**2)
e1 = c1/a1
e2 = c2/(a1**2 + a2)


def piCal(genotype_arr, num_sites):
    '''calculate the pi from parsed genotype arrays
    Note: assume the window size is 1.
    '''
    # preset variable
    num_sequences = len(genotype_arr)
    print(f'The number of sequences is {num_sequences}')

    pi = 0
    for i in range(0, len(genotype_arr)):
        for j in range(i+1, len(genotype_arr)):
            num_diffnul = sum ( genotype_arr[i][base]!=genotype_arr[j][base] for base in range(num_sites) )
            # print(num_diffnul)
            pi += num_diffnul
        # print(pi)
    pi /= (num_sequences*(num_sequences-1)*0.5)
    return pi
pi = piCal(genotype_arr, num_sites)
print(f'pi: {pi}')


def thetahatCal(segregating_site, a1):
    return segregating_site/a1
theta_hat = thetahatCal(segregating_site, a1)
print(f'theta_hat: {theta_hat}')


def tajimaD_Cal(pi, segregating_site, theta_hat, e1, e2):
    return (pi - theta_hat) / np.sqrt(e1*segregating_site + e2*segregating_site*(segregating_site-1))
D = tajimaD_Cal(pi, segregating_site, theta_hat, e1, e2)
print(f'Tajima\'s D: {D}')