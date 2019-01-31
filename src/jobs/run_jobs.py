import os

for file_name in os.listdir('.'):
	if file_name.endswith('sh'):
		cmd = 'sbatch --mem=12gb --qos throughput --time=10:00:00 --nodes=1 --ntasks=4 ' + file_name
		os.system(cmd)