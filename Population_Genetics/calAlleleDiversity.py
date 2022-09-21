'''this module is used to calculate statistics in population genetics'''
from itertools import product
from collections import Counter

# class VCFLine():
#     def __init__(self, alleles, genotypes):
#         self.alleles = alleles
#         self.genotypes = genotypes


def CalAllelicdDiversity(filename):
    '''this function is used to parse the vcf file and deposit info'''
    line_count = 0
    allele_counts = 0
    with open(filename, 'rt') as input:
        for line in input:
            if "#" in line: 
                continue
            ref = line.strip().split('\t')[3]
            alt = line.strip().split('\t')[4]

            # Parse REF and ALT
            new_alt = []
            if len(ref) >= 2:   # skip the line with no ref allele
                continue
            for index, i in enumerate(alt.split(','), start=1):
                # print(index)
                if len(i) == 1:
                    new_alt.append(index)

            genotype = [i for i in product([0], new_alt)]   # it should be saved in the list format or it will be released
            if len(genotype) == 0:
                continue
            # print(genotype)

            # Count corresponding allele types
            joined_genotype = []
            for i in genotype:
                joined_genotype.append(f"{str(i[0])}|{str(i[1])}")
                joined_genotype.append(f"{str(i[1])}|{str(i[0])}")
            # print(joined_genotype)

            sample_genotype = [i.split(':')[0] for i in line.strip().split('\t')[9:]]
            # print(sample_genotype)
            raw_allele_counts = Counter(sample_genotype)
            # print(raw_allele_counts)
            for key in list(raw_allele_counts.keys()):
                if key not in joined_genotype:
                    raw_allele_counts.pop(key, 'None')
            # print(raw_allele_counts)
            line_count += 1
            allele_counts += sum([i[1] for i in raw_allele_counts.items()])
            # print(allele_counts)
    return allele_counts*2 / line_count
                    

ad = CalAllelicdDiversity('pinfsc50_filtered.vcf')
print(ad)
# print(genotypes)

