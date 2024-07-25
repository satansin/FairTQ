'''
Created on Jan 24, 2017

@author: meike.zehlike
'''

import pandas as pd
from datasetCreator.candidate import Candidate

from time import process_time

from datasetCreator.preprocess import preprocessData
import os.path


'''
def createGender(filename, *columnsToRead):
    """
    currently working with recidivism score as qualification attribute in candidates. Change index
    to try with other columns
    """
    nonProtected = []
    protected = []
    data = None
    with open(filename) as csvfile:
        data = pd.read_csv(csvfile, usecols=columnsToRead)
        for row in data.itertuples():
            # change to different index in row[.] to access other columns from csv file
            if row[4] == 0:
                # nonProtected.append(Candidate(1 - row[3], [], id))
                nonProtected.append(Candidate(row[0] + row[1] + row[2], [], id))
            else:
                # protected.append(Candidate(1 - row[3], ["female"], id))
                protected.append(Candidate(row[0] + row[1] + row[2], ["female"], id))

    # sort candidates by recidivism scores in COMPAS
    protected.sort(key=lambda candidate: candidate.qualification, reverse=True)
    nonProtected.sort(key=lambda candidate: candidate.qualification, reverse=True)

    return protected, nonProtected, data
'''


def createRace(filename, dim):
    """
    currently working with recidivism score as qualification attribute in candidates. Change index
    to try with other columns
    """
    filepath = "../data_fairtq_processed/" + filename + ".csv"
    preprocessData(filename)
    # if not os.path.isfile(filepath):
    #     preprocessData(filename)

    cols = []
    for i in range(dim):
        cols.append("C" + str(i + 1))
    cols.append("P")

    nonProtected = []
    protected = []
    data = None
    with open(filepath) as csvfile:
        data = pd.read_csv(csvfile, usecols=cols)
        id = 1
        for row in data.itertuples():
            # print(row)
            # print(row[0], row[1], row[2], row[3])
            ### Note that row[0] is the index
            # change to different index in row[.] to access other columns from csv file
            qualification = 0
            for i in range(dim):
                qualification += row[i + 1]
            if row[dim + 1] == 0:
                nonProtected.append(Candidate(qualification, ["black"], id))
            else:
                ## note that we set non-black as protected, because we in general the black are discriminated if they are more selected in top-k
                protected.append(Candidate(qualification, [], id))
            id += 1

    t1_start = process_time()
    # sort candidates by decile scores in COMPAS
    protected.sort(key=lambda candidate: candidate.qualification, reverse=True)
    nonProtected.sort(key=lambda candidate: candidate.qualification, reverse=True)
    t1_stop = process_time()

    list_data = data.values.tolist()
    print("number of protected:", len(protected))
    print("example 1:", list_data[protected[0].id - 1])
    print("example 2:", list_data[protected[len(protected) - 1].id - 1])
    print("number of nonProtected:", len(nonProtected))
    print("example 1:", list_data[nonProtected[0].id - 1])
    print("example 2:", list_data[nonProtected[len(nonProtected) - 1].id - 1])

    return protected, nonProtected, data, (t1_stop - t1_start)
