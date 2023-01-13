'''this module is used to translate the sorted CDS sequence into protein sequence
Note: this module now cannot tackle partial codon problems
# ---------------------------------------------------------------------------
BiopythonWarning: Partial codon, len(sequence) not a multiple of three. 
Explicitly trim the sequence or add trailing N before translation. 
This may become an error in future.
# ---------------------------------------------------------------------------
'''
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def DNATranslate(inputname, outputname):
    '''
    Do DNA sequence translation all at once
    '''
    m_records = []
    # PROTEIN = open(output, 'w')
    with open(inputname, 'rt') as input:
        for record in SeqIO.parse(input, 'fasta'):
            header = record.description
            seq = record.seq.translate()
            
            # Direct write to output
            # PROTEIN.write(f'>{header}\n')
            # PROTEIN.write(f'{seq}\n')
            
            # Save as SeqRecords
            m_records.append(SeqRecord(
                Seq(seq),
                id=header,
                description=''
                ))
    
    SeqIO.write(m_records, outputname, "fasta")

def main():
    import optparse

    usage = """python -m plastidUtilis.Translate -i <sorted_CDS> -o <output_protein_sequence>
                                                                                -- Youpu=Chen"""
    parser = optparse.OptionParser(usage)
    parser.add_option("-i", dest="input", help="input sorted CDS fasta",
                    metavar="FILE", action="store", type="string")
    parser.add_option("-o", dest="output", help="output filename",
                    metavar="OUT", action="store", type="string")
    (options, args) = parser.parse_args()

    DNATranslate(options.input, options.output)


if __name__ == "__main__":
    main()


