'''this module is used to detect the abnormal CDS sequence annotated by Geseq
Note: input sequence should be in protein / amino acid format. 
'''
def AbnormalDetect(input, output):
    aa_dict = {}
    # aa_file = options.aa
    aa_file = input
    # out_file = options.output
    out_file = output
    with open(aa_file) as f:
        data = f.readlines()
        for line in data:
            if line.startswith(">"):
                aa_header = line.strip()
                aa_dict[aa_header] = []
            else:
                aa_dict[aa_header].append(line.strip())

    for k, v in aa_dict.items():
        aa_dict[k] = "".join(v)

    normal_dict = {}
    abnormal_dict = {}

    for k, v in aa_dict.items():
        # normal
        if v.endswith("*") and int(len(v.split("*"))) == 2:
            normal_dict[k.split()[0][1:]] = v
        else:
            # abnormal
            abnormal_dict[k.split()[0][1:]] = v
    with open(out_file, "w") as out:
        for k, v in abnormal_dict.items():
            if k in normal_dict.keys():
                # copy abnormal
                out.write(k+ "\t" + "copy_abnormal" + "\n")
            elif k not in normal_dict.keys():
                out.write(k + "\n")
    out.close()


def main():
    import optparse

    usage = """python -m plastidUtilis.AbnormalDetect -a <input_fasta> -o <output_filename>
                                                                                        --Joe"""
    parser = optparse.OptionParser(usage)
    parser.add_option("-a", dest="aa", help="input protein fasta",
                    metavar="FILE", action="store", type="string")
    parser.add_option("-o", dest="output", help="output filename",
                    metavar="OUT", action="store", type="string")
    (options, args) = parser.parse_args()

    AbnormalDetect(options.aa, options.output)


if __name__ == "__main__":
    main()
