import argparse
import json

parser = argparse.ArgumentParser()
parser.add_argument("result", help="The result file to clean")
args = parser.parse_args()

rot_table_1 = json.load(open('dna/table.json'))
rot_table_2 = json.load(open(args.result))

for key in rot_table_1:
    for i in range(3):
        inf = rot_table_1[key][i] - rot_table_1[key][3+i]
        sup = rot_table_1[key][i] + rot_table_1[key][3+i]
        if rot_table_2[key][i] < inf or rot_table_2[key][i] > sup:
            raise ValueError(f"Value {rot_table_2[key][i]} for {key} is out of range [{inf}, {sup}]")
        if (rot_table_2[key][3+i][0] != rot_table_2[key][i] - inf or rot_table_2[key][3+i][1] != sup - rot_table_2[key][i]):
            print('Correcting', key)
        rot_table_2[key][3+i][0] = rot_table_2[key][i] - inf
        rot_table_2[key][3+i][1] = sup - rot_table_2[key][i]

new_name = args.result.split('/')
new_name[-1] = 'cleaned_' + new_name[-1]
new_name = '/'.join(new_name)
with open(new_name, "w") as file:
    json.dump(rot_table_2, file, indent=4)