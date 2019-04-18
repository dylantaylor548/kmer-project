import os
import sys

files = []
with open(sys.argv[1]) as f:
	for line in f:
		val = line.strip()
		files.append(val)


for i in range(0, len(files)):
	# print (files)
	# print(i, files[i])
	files_copy = [elem for elem in files]
	del files_copy[i]
	# print (files_copy)
	print (files[i]+'\t'+','.join(files_copy))
	# break

