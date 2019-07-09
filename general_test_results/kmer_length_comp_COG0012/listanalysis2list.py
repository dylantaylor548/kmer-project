import os

for file_name in os.listdir(os.getcwd()):
	if file_name.endswith('.list_analysis'):
		fo = open(file_name[:-9],'w')
		with open(file_name) as opened_file:
			for line in opened_file:
				line = line.strip()
				if len(line.split(',')) == 2 and not line.startswith('runtime'):
					fo.write(line.split(',')[0] + '\n')
		fo.close()