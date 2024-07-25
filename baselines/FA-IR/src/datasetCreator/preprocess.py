import csv

def preprocessData(filename):
	filepath = "../data_fairtq/" + filename + ".txt"
	dim = 0
	data = []
	size = 0
	with open(filepath) as file:
		for line in file:
			# print(line.rstrip())
			items = line.rstrip().split()

			dim = (len(items) - 1) // 2

			line_data = []
			for i in range(dim):
				line_data.append((float(items[i]) + float(items[i + dim])) / 2)
			line_data.append(items[len(items) - 1])

			data.append(line_data)
			size = size + 1

	csvheader = []
	for i in range(dim):
		csvheader.append("C" + str(i + 1))
	csvheader.append("P")

	csvpath = "../data_fairtq_processed/" + filename + ".csv"
	with open(csvpath, 'w', newline='') as csvfile:
		writer = csv.writer(csvfile)
		writer.writerow(csvheader)
		for data_row in data:
			writer.writerow(data_row)