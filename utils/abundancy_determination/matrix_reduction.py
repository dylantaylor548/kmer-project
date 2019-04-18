from copy import deepcopy

def reduce_matrix(matrix):
	nrows = len(matrix)
	ncols = len(matrix[0])
	row_num = 0
	col_num = 0

	while (row_num < nrows) and (col_num < ncols): 

		test = False

		for i in range(row_num,nrows):
			row = deepcopy(matrix[i])
			if row[col_num] != 0:
				row = [int(x / row[col_num]) for x in row]
				hold = deepcopy(matrix[row_num])
				matrix[row_num] = row
				matrix[i] = hold
				test = True
				break

		for i in range(0,nrows):
			if i != row_num:
				if matrix[i][col_num] != 0:
					mult = matrix[i][col_num]
					mult_row = [(x*mult) for x in matrix[row_num]]
					matrix[i] = [matrix[i][x] - mult_row[x] for x in range(0,ncols)]


		col_num += 1
		if test == True:
			row_num += 1

	return matrix









matrix = [
[1,1,0,0,0,0,0,0,0,0,0,0,78],
[0,0,1,1,0,0,0,0,0,0,0,0,259],
[0,0,0,0,1,1,0,0,0,0,0,0,609],
[0,0,0,0,0,0,1,1,0,0,0,0,520],
[0,0,0,0,0,0,0,0,1,1,1,1,196],
[0,0,1,0,0,0,0,0,1,0,0,0,280],
[1,0,0,1,1,0,0,0,0,1,0,0,500],
[0,0,0,0,0,1,1,0,0,0,1,0,496],
[0,1,0,0,0,0,0,1,0,0,0,1,386]
]

red_matrix = reduce_matrix(matrix)

for row in red_matrix:
	print(row)
