'''this module is used to parse the Hierarchical_Orthogroups.Count.csv'''
import sys

def ParseCount(filename):
    group_dict = {}
    with open(filename, 'rt') as input:
        next(input) # skip header line
        for line in input:
            parsed_line = line.strip().split(',')
            for x in range(1, len(parsed_line)):
                
                if parsed_line[x] == 'NA':     # turn NA into a flag, e.g. -1 in string format
                    parsed_line[x] = '-1'

                elif int(parsed_line[x]) > 0:  # turn values > 0 into a flag, e.g. C
                    parsed_line[x] = 'C'

            # print(','.join(parsed_line[1:]))
            key = ','.join(parsed_line[1:])
            group_dict[key] = group_dict.get(key, 0) + 1
    return group_dict



group_dict = ParseCount(sys.argv[1])
# print(group_dict)
# print(len(group_dict.keys()))
output = open(sys.argv[2], 'w')
for i in group_dict.items():
    output.write(f'{i[0]}\t{i[1]}\n')
output.close()
