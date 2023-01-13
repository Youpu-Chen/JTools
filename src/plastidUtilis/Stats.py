'''this module is used to conduct basic statistics based on multiple bio-formats'''
from Bio import SeqIO

def SeqLength(filename):
    '''Calculate sequence length in a given fasta file'''
    with open(filename) as input:
        for record in SeqIO.parse(input, 'fasta'):
            print(len(record.seq))

def main():
    import optparse

    usage = """python -m plastidUtilis.Stats -i <sorted_CDS>
                                                            -- Youpu=Chen"""
    parser = optparse.OptionParser(usage)
    parser.add_option("-i", dest="input", help="input fasta",
                    metavar="FILE", action="store", type="string")
    (options, args) = parser.parse_args()

    SeqLength(options.input)


if __name__ == "__main__":
    main()