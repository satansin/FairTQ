import json
import os
import csv
import pdb

def get_num_months(dates):
	try:
		if len(dates) == 11:
			start_ye = int(dates[0:4])
			end_ye = int(dates[7:11])
			return (12 * (end_ye - start_ye + 1))
		elif len(dates) == 17:
			start_mo = int(dates[0:2])
			start_ye = int(dates[3:7])
			end_mo = int(dates[10:12])
			end_ye = int(dates[13:17])
			return (13 - start_mo + end_mo + 12 * (end_ye - start_ye - 1))
		elif len(dates) == 14:
			if dates[2:3] == "/":
				start_mo = int(dates[0:2])
				start_ye = int(dates[3:7])
				end_mo = 12
				end_ye = int(dates[10:14])
			else:
				start_mo = 1
				start_ye = int(dates[0:4])
				end_mo = int(dates[7:9])
				end_ye = int(dates[10:14])
			return (13 - start_mo + end_mo + 12 * (end_ye - start_ye - 1))
		else:
			# print(len(dates))
			return 0
	except ValueError:
		return 0


num_files = 57

raw_info = [["num_mo_job", "num_mo_edu", "sex"]]

max_num_mo_job = -1
max_num_mo_edu = -1

global_i = 0

for i in range(1, num_files + 1):
	start = (i - 1) * 40 + 1
	end = i * 40
	file_name = f"SHAano{i:02d}_{start}-{end}.json"

	with open(os.path.join(".", "Xing", file_name), encoding="utf8") as json_file:
		json_content = json.loads(json_file.read())

		# print(file_name)

		for j in range(len(json_content["profiles"])):
			# if "profile" in json_content["profiles"][j].keys():
			# 	print(f"len(json_content['profiles'][{j}]['profile']):", len(json_content["profiles"][j]["profile"]))
			# else:
			# 	print(f"len(json_content['profiles'][{j}]['profile']):", 0)
			# if "education" in json_content["profiles"][j].keys():
			# 	print(f"len(json_content['profiles'][{j}]['education']):", len(json_content["profiles"][j]["education"]))
			# else:
			# 	print(f"len(json_content['profiles'][{j}]['education']):", 0)

			global_i += 1

			sex = json_content["profiles"][j]["profile"][0]["sex"]

			num_mo_job = 0
			if "jobs" in json_content["profiles"][j]["profile"][0]:
				for job in json_content["profiles"][j]["profile"][0]["jobs"]:
					if "jobDates" in job.keys():
						num_mo_job += get_num_months(job["jobDates"])
						# print(job["jobDates"])

			if  max_num_mo_job < num_mo_job:
				max_num_mo_job = num_mo_job

			num_mo_edu = 0
			if "education" in json_content["profiles"][j]:
				for edu_item in json_content["profiles"][j]["education"]:
					if "eduDuration" in edu_item.keys():
						num_mo_edu += get_num_months(edu_item["eduDuration"])
						# if global_i == 268:
						# print(edu_item["eduDuration"])

			if  max_num_mo_edu < num_mo_edu:
				max_num_mo_edu = num_mo_edu

			raw_info.append([num_mo_job, num_mo_edu, sex])

			# if global_i == 268:
			# 	print("sex:", sex)
			# 	print("num_mo_job:", num_mo_job)
			# 	print("num_mo_edu:", num_mo_edu)

with open('xing_raw.csv', 'w', newline='') as raw_file:
	csvwriter = csv.writer(raw_file)
	csvwriter.writerows(raw_info)

processed_info = []
for i, raw_item in enumerate(raw_info):
	if i == 0:
		continue

	num_mo_job_norm = round(raw_item[0] / max_num_mo_job, 4)
	num_mo_edu_norm = round(raw_item[1] / max_num_mo_edu, 4)
	sex_norm = 1 if raw_item[2] == "f" else 0
	processed_info.append([num_mo_job_norm + 0.0001, num_mo_edu_norm + 0.0001, num_mo_job_norm - 0.0001, num_mo_edu_norm - 0.0001, sex_norm])


with open('xing_2d_2280.txt', 'w', newline='') as raw_file:
	csvwriter = csv.writer(raw_file, delimiter=' ',)
	csvwriter.writerows(map(lambda t: ("%.4f" % t[0], "%.4f" % t[1], "%.4f" % t[2], "%.4f" % t[3], t[4]), processed_info))