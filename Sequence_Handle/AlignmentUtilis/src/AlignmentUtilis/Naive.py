import random
import subprocess
import os

class Naive():
    @staticmethod
    def read_Cicurlar_Genome(filename):
        '''this function is used to read sequence file in fasta format
        Note: it could only be applied to read one sequence, which could be the genome of chloroplast, mitochondia as well as microbes'''
        genome = ''
        with open(filename, 'r') as f:
            for line in f:
                # ignore header line with genome information
                if not line[0] == '>':
                    genome += line.rstrip()
        return genome
    
    @staticmethod
    def readFastq(filename):
        '''this function is used to read fastq file (not in gzip or any compressed format), 
        it return sequence and its correspondng sequence quality'''
        sequence = []
        qualities = []
        with open(filename) as fh:
            while True:
                # in a loop, it will run codes below, and read only a line each time
                fh.readline()
                seq = fh.readline().rstrip()
                fh.readline()
                qual = fh.readline().rstrip()

                if len(seq) == 0:
                    break
                sequence.append(seq)
                qualities.append(qual)
        return sequence, qualities

    @staticmethod
    def reverseComplement(s):
        complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'} # N represent the ambiguous bases
        t = ''
        for base in s:
            t = complement[base] + t   # pre-pending -> reverse the string
        return t
    
    @staticmethod
    def naive_exact_matching(pattern, target):
        '''this function is the application of naive exact matching algorithms to count the occurence of input sequence
        Note: this version has not included mismatches'''
        occurrences = []
        for i in range(len(target) - len(pattern) + 1):
            match = True
            for j in range(len(pattern)):
                if target[i+j] != pattern[j]:
                    match = False
                    break
            if match:
                occurrences.append(i)
        return occurrences
    @staticmethod
    def generateReads(genome, numReads, readLen):
        ''' this function return a set of reads from a given genome'''
        reads = []
        for _ in range(numReads):
            start = random.randint(0, len(genome)-readLen) - 1
            reads.append(genome[start : start+readLen])
        return reads

    def __init__(self, reads, genome):
        self.reads = reads
        self.genome = genome  # The reverse complement should be considered, as for the most reliable results

    def runNaive(self, seed_flag=True, reverse_flag=True, seed_length=30):
        '''this function is the main function to run naive exact matching
        Note: it's recommanded to run the alignment process against the reverse complement of genome'''
        numMatched = 0
        n = 0
        # Using custom seed length and run reverse complement alignment
        if seed_flag == True and reverse_flag == True:
            for read in self.reads:
                read = read[:seed_length]
                # print(read)
                matches = Naive.naive_exact_matching(read, self.genome)
                matches.extend(Naive.naive_exact_matching(Naive.reverseComplement(read), self.genome))
                n += 1
                if len(matches) > 0:
                    numMatched += 1
            print('%d / %d reads matched the genome exactly!' % (numMatched, n))

        # Using whole sequence length and run reverse complement alignment
        if seed_flag == False and reverse_flag == True:
            for read in self.reads:
                # read = read[:seed_length]
                matches = Naive.naive_exact_matching(read, self.genome)
                matches.extend(Naive.naive_exact_matching(Naive.reverseComplement(read), self.genome))
                n += 1
                if len(matches) > 0:
                    numMatched += 1
            print('%d / %d reads matched the genome exactly!' % (numMatched, n))


def runCMD(cmd, verbose = False, *args, **kwargs):
    '''this function is used to run commandlines e.g. wget, curl.
    Note: from Roel Peters, "Using Python and wget to Download Web Pages and Files" '''
    process = subprocess.Popen(
        cmd,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE,
        text = True, 
        shell=True
    )
    std_out, std_err = process.communicate()
    if verbose:
        print(std_out.strip(), std_err)
    pass

def getExampleDatasets():
    '''this function is used to helpe user to retrieve example datasets.
    Note: only executed when the test dataset does not exist'''
    runCMD('mkdir exampledata')
    runCMD('wget http://d28rh4a8wq0iu5.cloudfront.net/ads1/data/phix.fa -O exampledata/phix.fa')
    runCMD('wget http://d28rh4a8wq0iu5.cloudfront.net/ads1/data/ERR266411_1.first1000.fastq -O exampledata/ERR266411_1.first1000.fastq')


if __name__ == "__main__":
    # Basic utility test
    print('------------------------------- Naive Exact Macthing Basic Utility Test -------------------------------')
    test_occurrences = Naive.naive_exact_matching('AG', 'AGCTTAGATAGC')
    print('The pattern is AG')
    print('The target sequence is AGCTTAGATAGC')
    print(f'The start position of exact matching is {test_occurrences}')
    
    # Integrated utility test
    getExampleDatasets()   # automatically download the datasets

    # 1. Using actual reads
    print('------------------------------- Naive Exact Macthing: Using Actual Reads -------------------------------')
    phix_genome = Naive.read_Cicurlar_Genome('exampledata/phix.fa')
    phix_reads, _ = Naive.readFastq('exampledata/ERR266411_1.first1000.fastq')
    phix_naive_real = Naive(phix_reads, phix_genome)
    phix_naive_real.runNaive(seed_flag=True, reverse_flag=True, seed_length=30)

    # 2. Using simulated reads
    print('------------------------------- Naive Exact Macthing Basic Utility Test: Using Simulated Reads -------------------------------')
    simulated_reads = Naive.generateReads(phix_genome, 100, 100)  # generate 100 reads in 100 bp
    phix_naive_sim = Naive(simulated_reads, phix_genome)
    phix_naive_sim.runNaive(seed_flag=True, reverse_flag=True, seed_length=30)