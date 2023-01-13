'''this module is parse the self-made TABLE file into GenBank format
Example TABLE: it should contain features like: 
1) nad1, nad2, nad5 gene: trans-splicing gene
#nad1
nad1_1	gene	361447	361833	-
nad1_1	CDS	361447	361833	-	1
nad1_1	note	trans_splicing	trans_splicing

2) tRNA gene
3) rRNA gene

!!! All the product should be filled before the GBK transformation
Dev: What I want to do is automatically filling!!!
'''
def TABLE2GBK(input, output):
    gene_dict = {}
    # with open(options.table) as f:
    with open(input, 'rt', encoding='utf-8') as f:   
        # if not set "utf-8", it will report 
        # UnicodeDecodeError: 'gbk' codec can't decode byte 0xaa in position 10: illegal multibyte sequence

        data = f.readlines()
        gene_num = 0
        for line in data:
            # not blank
            if line != "\n" and not line.startswith("#"):
                treat = line.split("\t")[1]
                if treat == 'gene':
                    gene_name = line.split("\t")[0] + "@" + str(gene_num)
                    gene_dict[gene_name] = {}
                    # start end orientation
                    gene_dict[gene_name]['gene'] = [line.split()[2], line.split()[3], line.split()[4]]
                    gene_num = gene_num + 1
                    continue
                elif treat == 'CDS' or treat == 'tRNA' or treat == 'rRNA':
                    if treat not in gene_dict[gene_name].keys():
                        gene_dict[gene_name][treat] = {}
                        exon = line.split()[5]
                        # start end orientation
                        gene_dict[gene_name][treat][exon] = [line.split()[2], line.split()[3], line.split()[4]]
                        continue
                    elif treat in gene_dict[gene_name].keys():
                        exon = line.split()[5]
                        # start end orientation
                        gene_dict[gene_name][treat][exon] = [line.split()[2], line.split()[3], line.split()[4]]
                        continue

                elif treat == 'note':
                    note = line.split("\t")[2]
                    print(note)
                    content = line.split("\t")[3].strip()
                    gene_dict[gene_name][note] = content
            

    # print gbk
    front_gap = " " * 5
    # with open(out_gbk, "w") as out:
    with open(output, "w") as out:
        for k, v in gene_dict.items():
            # gene
            rear_gap = " " * (21 - 5 - int(len("gene")))
            start = gene_dict[k]["gene"][0]
            end = gene_dict[k]["gene"][1]
            sign = gene_dict[k]["gene"][2]
            out.write(front_gap + "gene" + rear_gap)
            if sign == "+":
                out.write(start + ".." + end + "\n")
                out.write(21 * " " + '/gene="%s"' % k.split("@")[0] + "\n")
            elif sign == "-":
                out.write("complement(" + start + ".." + end + ")" + "\n")
                out.write(21 * " " + '/gene="%s"' % k.split("@")[0] + "\n")

            # CDS
            if "CDS" in gene_dict[k].keys() and "trans_splicing" not in gene_dict[k].keys():
                rear_gap = " " * (21 - 5 - int(len("CDS")))
                out.write(front_gap + "CDS" + rear_gap)
                # no intron
                exon_num = int(len(gene_dict[k]["CDS"].keys()))
                if exon_num == 1:
                    if sign == "+":
                        out.write(start + ".." + end + "\n")
                        out.write(21 * " " + '/gene="%s"' % k.split("@")[0] + "\n")
                    elif sign == "-":
                        out.write("complement(" + start + ".." + end + ")" + "\n")
                        out.write(21 * " " + '/gene="%s"' % k.split("@")[0] + "\n")
                # intron
                elif exon_num > 1:
                    exon_list = []
                    if sign == "+":
                        for i in range(1, exon_num + 1, 1):
                            start = gene_dict[k]["CDS"][str(i)][0]
                            end = gene_dict[k]["CDS"][str(i)][1]
                            exon_list.append(start + ".." + end)
                        out.write("join(" + ",".join(exon_list) + ")" + "\n")
                    elif sign == "-":
                        for i in range(exon_num, 0, -1):
                            start = gene_dict[k]["CDS"][str(i)][0]
                            end = gene_dict[k]["CDS"][str(i)][1]
                            exon_list.append(start + ".." + end)
                        out.write("complement(join(" + ",".join(exon_list) + "))" + "\n")
                    out.write(21 * " " + '/gene="%s"' % k.split("@")[0] + "\n")

                # product
                if "product" in gene_dict[k].keys():
                    out.write(21 * " " + '/product="%s"' % gene_dict[k]["product"] + "\n")

                # translation
                if "translation" in gene_dict[k].keys():
                    translation = '/translation="%s"' % gene_dict[k]["translation"]
                    translation_list = [translation[i: i + 58] for i in range(0, len(translation), 58)]
                    for i in range(0, len(translation_list)):
                        out.write(21 * " " + translation_list[i] + "\n")

                # exon and intron
                exon_rear_gap = " " * (21 - 5 - int(len("exon")))
                intron_rear_gap = " " * (21 - 5 - int(len("intron")))
                if exon_num > 1:
                    if sign == "+":
                        for i in range(1, exon_num + 1, 1):
                            start = gene_dict[k]["CDS"][str(i)][0]
                            end = gene_dict[k]["CDS"][str(i)][1]
                            out.write(front_gap + "exon" + exon_rear_gap + start + ".." + end + "\n")
                            out.write(21 * " " + '/gene="%s"' % k.split("@")[0] + "\n")
                            out.write(21 * " " + '/number="%s"' % str(i) + "\n")
                            intron_num = i
                            if intron_num < exon_num:
                                intron_start = int(end) + 1
                                intron_end = int(gene_dict[k]["CDS"][str(i + 1)][0]) - 1
                                out.write(
                                    front_gap + "intron" + intron_rear_gap + str(intron_start) + ".." + str(
                                        intron_end) + "\n")
                                out.write(21 * " " + '/gene="%s"' % k.split("@")[0] + "\n")
                                out.write(21 * " " + '/number="%s"' % str(intron_num) + "\n")

                    if sign == "-":
                        for i in range(exon_num, 0, -1):
                            start = gene_dict[k]["CDS"][str(i)][0]
                            end = gene_dict[k]["CDS"][str(i)][1]
                            out.write(front_gap + "exon" + exon_rear_gap + "complement(" + start + ".." + end + ")" + "\n")
                            out.write(21 * " " + '/gene="%s"' % k.split("@")[0] + "\n")
                            out.write(21 * " " + '/number="%s"' % str(i) + "\n")
                            intron_num = i - 1
                            if intron_num > 0:
                                intron_start = int(end) + 1
                                intron_end = int(gene_dict[k]["CDS"][str(i - 1)][0]) - 1
                                out.write(
                                    front_gap + "intron" + intron_rear_gap + "complement(" + str(intron_start) + ".." + str(
                                        intron_end) + ")" + "\n")
                                out.write(21 * " " + '/gene="%s"' % k.split("@")[0] + "\n")
                                out.write(21 * " " + '/number="%s"' % str(intron_num) + "\n")

            # trans_splicing
            elif "CDS" in gene_dict[k].keys() and "trans_splicing" in gene_dict[k].keys():
                exon_rear_gap = " " * (21 - 5 - int(len("exon")))
                intron_rear_gap = " " * (21 - 5 - int(len("intron")))
                exon_num = int(len(gene_dict[k]["CDS"].keys()))
                # one exon
                if exon_num == 1:
                    out.write(front_gap + "exon" + exon_rear_gap)
                    if sign == "+":
                        out.write(start + ".." + end + "\n")
                    elif sign == "-":
                        out.write("complement(" + start + ".." + end + ")" + "\n")
                    out.write(21 * " " + '/gene="%s"' % k.split("@")[0] + "\n")
                    out.write(21 * " " + '/number="%s"' % str(list(gene_dict[k]["CDS"].keys())[0]) + "\n")

                # exon and intron
                if exon_num > 1:
                    intron_num = 1
                    if sign == "+":
                        num_list = list(sorted(list([gene_dict[k]["CDS"].keys()])[0]))
                        for i in num_list:
                            start = gene_dict[k]["CDS"][str(i)][0]
                            end = gene_dict[k]["CDS"][str(i)][1]
                            out.write(front_gap + "exon" + exon_rear_gap + start + ".." + end + "\n")
                            out.write(21 * " " + '/gene="%s"' % k.split("@")[0] + "\n")
                            out.write(21 * " " + '/number="%s"' % str(i) + "\n")
                            if intron_num < exon_num:
                                intron_start = int(end) + 1
                                intron_end = int(gene_dict[k]["CDS"][str(int(i) + 1)][0]) - 1
                                out.write(
                                    front_gap + "intron" + intron_rear_gap + str(intron_start) + ".." + str(
                                        intron_end) + "\n")
                                out.write(21 * " " + '/gene="%s"' % k.split("@")[0] + "\n")
                                out.write(21 * " " + '/number="%s"' % str(intron_num) + "\n")
                            intron_num = intron_num + 1

                    if sign == "-":
                        num_list = list(sorted(list([gene_dict[k]["CDS"].keys()])[0], reverse=True))
                        for i in num_list:
                            start = gene_dict[k]["CDS"][str(i)][0]
                            end = gene_dict[k]["CDS"][str(i)][1]
                            out.write(front_gap + "exon" + exon_rear_gap + "complement(" + start + ".." + end + ")" + "\n")
                            out.write(21 * " " + '/gene="%s"' % k.split("@")[0] + "\n")
                            out.write(21 * " " + '/number="%s"' % str(i) + "\n")
                            if intron_num < exon_num:
                                intron_start = int(end) + 1
                                intron_end = int(gene_dict[k]["CDS"][str(int(i) - 1)][0]) - 1
                                out.write(
                                    front_gap + "intron" + intron_rear_gap + "complement(" + str(intron_start) + ".." + str(
                                        intron_end) + ")" + "\n")
                                out.write(21 * " " + '/gene="%s"' % k.split("@")[0] + "\n")
                                out.write(21 * " " + '/number="%s"' % str(intron_num) + "\n")
                            intron_num = intron_num + 1


            # tRNA
            elif "tRNA" in gene_dict[k].keys():
                tRNA_rear_gap = " " * (21 - 5 - int(len("tRNA")))
                out.write(front_gap + "tRNA" + tRNA_rear_gap)
                # no intron
                exon_num = int(len(gene_dict[k]["tRNA"].keys()))
                if exon_num == 1:
                    if sign == "+":
                        out.write(start + ".." + end + "\n")
                        out.write(21 * " " + '/gene="%s"' % k.split("@")[0] + "\n")
                    elif sign == "-":
                        out.write("complement(" + start + ".." + end + ")" + "\n")
                        out.write(21 * " " + '/gene="%s"' % k.split("@")[0] + "\n")
                # intron
                elif exon_num > 1:
                    exon_list = []
                    if sign == "+":
                        for i in range(1, exon_num + 1, 1):
                            start = gene_dict[k]["tRNA"][str(i)][0]
                            end = gene_dict[k]["tRNA"][str(i)][1]
                            exon_list.append(start + ".." + end)
                        out.write("join(" + ",".join(exon_list) + ")" + "\n")
                    elif sign == "-":
                        for i in range(exon_num, 0, -1):
                            start = gene_dict[k]["tRNA"][str(i)][0]
                            end = gene_dict[k]["tRNA"][str(i)][1]
                            exon_list.append(start + ".." + end)
                        out.write("complement(join(" + ",".join(exon_list) + "))" + "\n")
                    out.write(21 * " " + '/gene="%s"' % k.split("@")[0] + "\n")
                    
                # product
                if "product" in gene_dict[k].keys():
                    out.write(21 * " " + '/product="%s"' % gene_dict[k]["product"] + "\n")

                # exon and intron
                exon_rear_gap = " " * (21 - 5 - int(len("exon")))
                intron_rear_gap = " " * (21 - 5 - int(len("intron")))
                if exon_num > 1:
                    if sign == "+":
                        for i in range(1, exon_num + 1, 1):
                            start = gene_dict[k]["tRNA"][str(i)][0]
                            end = gene_dict[k]["tRNA"][str(i)][1]
                            out.write(front_gap + "exon" + exon_rear_gap + start + ".." + end + "\n")
                            out.write(21 * " " + '/gene="%s"' % k.split("@")[0] + "\n")
                            out.write(21 * " " + '/number="%s"' % str(i) + "\n")
                            intron_num = i
                            if intron_num < exon_num:
                                intron_start = int(end) + 1
                                intron_end = int(gene_dict[k]["tRNA"][str(i + 1)][0]) - 1
                                out.write(
                                    front_gap + "intron" + intron_rear_gap + str(intron_start) + ".." + str(
                                        intron_end) + "\n")
                                out.write(21 * " " + '/gene="%s"' % k.split("@")[0] + "\n")
                                out.write(21 * " " + '/number="%s"' % str(intron_num) + "\n")

                    if sign == "-":
                        for i in range(exon_num, 0, -1):
                            start = gene_dict[k]["tRNA"][str(i)][0]
                            end = gene_dict[k]["tRNA"][str(i)][1]
                            out.write(front_gap + "exon" + exon_rear_gap + "complement(" + start + ".." + end + ")" + "\n")
                            out.write(21 * " " + '/gene="%s"' % k.split("@")[0] + "\n")
                            out.write(21 * " " + '/number="%s"' % str(i) + "\n")
                            intron_num = i - 1
                            if intron_num > 0:
                                intron_start = int(end) + 1
                                intron_end = int(gene_dict[k]["tRNA"][str(i - 1)][0]) - 1
                                out.write(
                                    front_gap + "intron" + intron_rear_gap + "complement(" + str(intron_start) + ".." + str(
                                        intron_end) + ")" + "\n")
                                out.write(21 * " " + '/gene="%s"' % k.split("@")[0] + "\n")
                                out.write(21 * " " + '/number="%s"' % str(intron_num) + "\n")


            # rRNA
            elif "rRNA" in gene_dict[k].keys():
                rRNA_rear_gap = " " * (21 - 5 - int(len("rRNA")))
                out.write(front_gap + "rRNA" + rRNA_rear_gap)
                if sign == "+":
                    out.write(start + ".." + end + "\n")
                    out.write(21 * " " + '/gene="%s"' % k.split("@")[0] + "\n")
                elif sign == "-":
                    out.write("complement(" + start + ".." + end + ")" + "\n")
                # product
                if "product" in gene_dict[k].keys():
                    out.write(21 * " " + '/product="%s"' % gene_dict[k]["product"] + "\n")
                    
            out.write("\n")
    out.close()


def main():
    import optparse
    usage = """python -m plastidUtilis.Table2GBK -t <self-made table file> -o <.gb> 
                                                                      --Joe"""

    parser = optparse.OptionParser(usage)
    parser.add_option("-t", dest="table", help="input annotation table",
                    metavar="FILE", action="store", type="string")
    parser.add_option("-o", dest="outname", help="output filename",
                    metavar="OUT", action="store", type="string")
    (options, args) = parser.parse_args()

    TABLE2GBK(options.table, options.outname)

if __name__ == "__main__":
    main()