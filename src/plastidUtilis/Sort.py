'''this module is used to do something with fasta
Dev: 
1) parse plastid genome and GFF3 file annotated from GeSeq and concatenate CDS sequences together into a whole-length gene sequence 
Only if the gene has multiple CDS, or it will just save one CDS sequence.
2) DNA sequence translation
'''
def ParseFasta(filename):
    '''Parse the fasta file with header , generated by TBtools.
    TBtools:
    java -cp <Path_to_TBtools>/TBtools_JRE1.6.jar biocjava.bioIO.GFF.ExtractFeaturefromGFF3andGenome \
    --inGtf Geseq.gff3 \
    --inGenome plastid.assembly \
    --outFile <output>.fasta \
    --targetFeature <selected_features, e.g. > --targetIdTag ID --retainAttr true --onlyCheck false --maxFeatureCounts 10 --minFeatureCounts 1 --upStreamBases 0 --downStreamBases 0 --onlyRetainFlank false --addN false
    '''
    # with open(options.Geseq_out, "r") as f:
    with open(filename, "r") as f:
        data = f.readlines()
        chr_gene_pos = {}
        gene_seq = {}
        chr_list = []
        for line in data:
            if line.startswith(">"):
                seq_header = line.strip()
                seq_info = line.split()
                for i in seq_info:
                    if "Location" in i:
                        pos = i
                        seqchr = pos.split("=")[1].split(":")[0].strip()
                        chr_list.append(seqchr)
                        seqstart = pos.split(":")[1]
                        chr_gene_pos.setdefault(seqchr, {})
                chr_gene_pos[seqchr][seq_header] = seqstart
                gene_seq[seq_header] = []
            else:
                gene_seq[seq_header].append(line.strip())
    return chr_gene_pos, gene_seq, chr_list

def SortByChromosome(chr_gene_pos, gene_seq, chr_list, outname, header_len):
    '''Sort the retrieved sequence based on their chr_number and start position'''
    chr_list = set(chr_list)
    if len(chr_list) > 1:
        sort_seq_chr = sorted(chr_list, key = lambda d:int(d[header_len:]), reverse = False)
    else:
        sort_seq_chr = chr_list
    # with open(options.outname, "w") as out:
    with open(outname, "w") as out:
        for chr in sort_seq_chr:
            # sort by start
            sort_seq_pos = sorted(chr_gene_pos[chr].items(), key=lambda d: int(d[1]), reverse=False)
            for i in sort_seq_pos:
                sort_seqname = i[0]
                seq = "".join(gene_seq[sort_seqname])
                strand = sort_seqname.split()[1].strip()
                for j in sort_seqname.split():
                    if "Len=" in j:
                        seqlen = j.split("=")[1]
                    if "Location" in j:
                        start = j.split(":")[1]
                        chr = j.split("=")[1].split(":")[0].strip()
                        end = int(int(start) + int(seqlen) - 1)
                    if "gene=" in j:
                        sim_seq_name = ">" + j.split("=")[1]
                out.write(sim_seq_name + "\t" + str(chr) + "\t" + str(start) + "\t" + str(end) + "\t" + strand + "\t" + str(seqlen) + "\n")
                out.write(seq + "\n")
    out.close()


def main():
    import optparse
    usage = """python -m plastidUtilis.Sort -f <Geseq_out_seq> -o <name_of_output> --header <Input_header>
                                                                            -Joe"""

    parser = optparse.OptionParser(usage)
    parser.add_option("-f", dest="Geseq_out", help="cds seqence extract from GeSeq output gff3",
                    metavar="FILE", action="store", type="string")
    parser.add_option("--header", dest="header", help="header of genome sequence, except the number",
                    metavar="STRING", action="store", type="string")
    parser.add_option("-o", dest="outname", help="output file name",
                    metavar="FILE", action="store", type="string")
    (options, args) = parser.parse_args()
    header_len = int(len(options.header))
    
    chr_gene_pos, gene_seq, chr_list = ParseFasta(options.Geseq_out)
    SortByChromosome(chr_gene_pos, gene_seq, chr_list, options.outname, header_len)


if __name__ == "__main__":
    main()