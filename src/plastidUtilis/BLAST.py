'''this module is able to execute local-installed BLAST to conduct sequence alignemnt analysis with preset parameters'''
import optparse
import os


def RunBLAST(INPUT, BLASTDATABASE, OUTPUT):
    trn_dict = {}
    # with open(tRNA_file) as f:
    with open(INPUT) as f:
        data = f.readlines()
        trn_dict["trash"] = []
        for line in data:
            if line.startswith(">"):
                seqname = line.strip().replace("\t", "@")[1:]
                if seqname not in trn_dict.keys():
                    trn_dict[seqname] = []
                elif seqname in trn_dict.keys():
                    seqname = "trash"
            else:
                trn_dict[seqname].append(line.strip())

    del trn_dict['trash']
    with open("temp_tRNA.fasta", "w") as t:
        for k, v in trn_dict.items():
            t.write(">" + k + "\n")
            t.write("".join(v) + "\n")
    t.close()

    print(len(trn_dict.keys()))
    os.system("makeblastdb -in %s -dbtype nucl -out %s" % (BLASTDATABASE, BLASTDATABASE))
    os.system(
        "blastn -query temp_tRNA.fasta -db %s -out blast.txt -perc_identity 90 -max_target_seqs 1 -subject_besthit  -num_threads 6 -outfmt  '6 qseqid qlen sseqid slen pident length score' " % database_file)
    os.system("rm temp_tRNA.fasta")

    blast_dict = {}
    with open("blast.txt") as b:
        data = b.readlines()
        for line in data:
            qseqid = line.split()[0]
            qlen = line.split()[1]
            sseqid = line.split()[2]
            slen = line.split()[3]
            pident = line.split()[4]
            length = line.split()[5]
            score = line.split()[6].strip()
            # length constrain
            if float(int(qlen) / int(slen)) > 0.9 and float(int(qlen) / int(slen)) < 1.1:
                # alignment length constrain
                if float(int(length) / int(slen)) > 0.9 and float(pident) > 0.90:
                    blast_dict[qseqid] = sseqid.split("-")[1] + "-" + sseqid.split("-")[2][:3]

    with open(OUTPUT, "w") as out:
        for k, v in blast_dict.items():
            old_id = k.split("@")[0]
            chr = k.split("@")[1]
            start = k.split("@")[2]
            end = k.split("@")[3]
            chain = k.split("@")[4]
            length = k.split("@")[5]
            new_id = v
            # intron
            if int(end) - int(start) + 1 == int(length):
                info = [new_id, chr, start, end, chain, length]
            else:
                info = [new_id, chr, start, end, chain, length, "intron"]
            out.write(">" + "\t".join(info) + "\n")
            out.write("".join(trn_dict[k]) + "\n")
    out.close()


def main():
    usage = """python -m plastidUtilis -d <tRNA_database> -t <sort_tRNA> -o <output.fasta>
                                                                                    --Joe"""
    parser = optparse.OptionParser(usage)
    parser.add_option("-d", dest="tRNA_database", help="tRNA blast database",
                    metavar="FILE", action="store", type="string")
    parser.add_option("-t", dest="sort_tRNA", help="sort_tRNA.fasta",
                    metavar="OUT", action="store", type="string")
    parser.add_option("-o", dest="output", help="output file name",
                    metavar="OUT", action="store", type="string")
    (options, args) = parser.parse_args()

    database_file = options.tRNA_database
    tRNA_file = options.sort_tRNA
    out_file = options.output

    RunBLAST(tRNA_file, database_file, out_file)


if __name__ == "__main__":
    main()