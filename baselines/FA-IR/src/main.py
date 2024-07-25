'''
Created on Mar 29, 2017

@author: meike.zehlike
'''
import os, sys
import argparse

import pdb

from time import process_time

from readWriteRankings.readAndWriteRankings import writePickleToDisk
from post_processing_methods.fair_ranker.create import fairRanking, feldmanRankingFunc
from utilsAndConstants.constants import ESSENTIALLY_ZERO
# from utilsAndConstants.utils import setMemoryLimit
from datasetCreator import compasData, xingData, allstarData, dotData, syntheticData #, germanCreditData, satData, xingProfilesReader, chileSATData
from evaluator.evaluator import Evaluator
from evaluator.failProbabilityYangStoyanovich import determineFailProbOfGroupFairnessTesterForStoyanovichRanking

EVALUATE_FAILURE_PROBABILITY = 0

def main():
    # setMemoryLimit(10000000000)

    # create the top-level parser
    parser = argparse.ArgumentParser(prog='FA*IR', description='a fair Top-k ranking algorithm',
                                     epilog="=== === === end === === ===")
    parser.add_argument("-c", "--create", nargs='*', help="creates a ranking from the raw data and dumps it to disk")
    parser.add_argument("-e", "--evaluate", nargs='*', help="evaluates rankings and writes results to disk")
    subparsers = parser.add_subparsers(help='sub-command help')

    # create the parser for the "create" command
    parser_create = subparsers.add_parser('dataset_create', help='choose a dataset to generate')
    parser_create.add_argument(dest='dataset_to_create', choices=["sat", "compas", "germancredit", "xing", "csat"])

    # create the parser for the "evaluate" command
    parser_evaluate = subparsers.add_parser('dataset_evaluate', help='choose a dataset to evaluate')
    parser_evaluate.add_argument(dest='dataset_to_evaluate', choices=["sat", "xing"
                                                                  "compas_gender", "compas_race",
                                                                  "germancredit_25", "germancredit_35", "germancredit_gender"])

    args = parser.parse_args()

    if args.create == []:
        print("creating rankings for all datasets...")
        createDataAndRankings()
    elif args.create == ['xing']:
        createAndRankXingData()
    elif args.create == ['compas']:
        createAndRankCOMPASData()
    elif args.create == ['dot']:
        createAndRankDOTData()
    elif args.create == ['allstar']:
        createAndRankAllstarData()
    elif args.create == ['synthetic']:
        createAndRankSyntheticData()
    # elif args.create == ['germancredit']:
    #     createAndRankGermanCreditData(100)
    # elif args.create == ['csat']:
    #     createAndRankChileData(1500)
    #=======================================================
    # elif args.evaluate == []:
    #     evaluator = Evaluator()
    #     evaluator.printResults()
    # elif args.evaluate == ['compas_gender']:
    #     evaluator = Evaluator('compas_gender')
    #     evaluator.printResults()
    #     evaluator.plotOrderingUtilityVsPercentageOfProtected()
    # elif args.evaluate == ['compas_race']:
    #     evaluator = Evaluator('compas_race')
    #     evaluator.printResults()
    # elif args.evaluate == ['germancredit_25']:
    #     evaluator = Evaluator('germancredit_25')
    #     evaluator.printResults()
    #     evaluator.plotOrderingUtilityVsPercentageOfProtected()
    # elif args.evaluate == ['germancredit_35']:
    #     evaluator = Evaluator('germancredit_35')
    #     evaluator.printResults()
    # elif args.evaluate == ['germancredit_gender']:
    #     evaluator = Evaluator('germancredit_gender')
    #     evaluator.printResults()
    # elif args.evaluate == ['xing']:
    #     evaluator = Evaluator('xing')
    #     evaluator.printResults()
    # elif args.evaluate == ['sat']:
    #     evaluator = Evaluator('sat')
    #     evaluator.printResults()

    else:
        print("FA*IR \n running the full program \n Press ctrl+c to abort \n \n")
        createDataAndRankings()
        # evaluator = Evaluator()
        # evaluator.printResults()

        # if EVALUATE_FAILURE_PROBABILITY:
        #     determineFailProbOfGroupFairnessTesterForStoyanovichRanking()


def createDataAndRankings():
    createAndRankXingData()
    createAndRankCOMPASData()
    createAndRankDOTData()
    createAndRankAllstarData()
    createAndRankSyntheticData()


# // std::string dataset;
# // int dataset_size;
# // int k;
# // int dim;
# // int num_groups;
# // float alpha;
# // float beta;
# // int one_div_theta;
# // int num_levels;
# // int rand_num_total;


def createAndRankXingData():
    pairsOfPAndAlpha = [(0.1, 0.1),  # no real results, skip in evaluation
                        (0.2, 0.1),  # no real results, skip in evaluation
                        (0.3, 0.1),  # no real results, skip in evaluation
                        (0.4, 0.1),  # no real results, skip in evaluation
                        (0.5, 0.0168),
                        (0.6, 0.0321),
                        (0.7, 0.0293),
                        (0.8, 0.0328),
                        (0.9, 0.0375)]

    exps_xing = [
        ("xing_2d_2280", 10, 2, 2),
        ("xing_2d_2280", 20, 2, 2),
        ("xing_2d_2280", 30, 2, 2),
        ("xing_2d_2280", 40, 2, 2),
        ("xing_2d_2280", 50, 2, 2),
        ("xing_2d_2280", 22, 2, 2),
    ]

    for exp_xing in exps_xing:
        protectedXingGender, nonProtectedXingGender, data, sort_time = xingData.createGender(exp_xing[0], exp_xing[2])
        rankAndDump(protectedXingGender, nonProtectedXingGender, exp_xing, pairsOfPAndAlpha, data, sort_time)


def createAndRankCOMPASData():

    pairsOfPAndAlpha = [(0.1, 0.0140),
                        (0.2, 0.0115),
                        (0.3, 0.0103),
                        (0.4, 0.0099),
                        (0.5, 0.0096),
                        (0.6, 0.0093),
                        (0.7, 0.0094),
                        (0.8, 0.0095),
                        (0.9, 0.0100)]

    exps_compas = [
        ("compas_3d_10k_2_subset_1000", 10, 3, 2),
        ("compas_3d_10k_2_subset_2000", 10, 3, 2),
        ("compas_3d_10k_2_subset_3000", 10, 3, 2),
        ("compas_3d_10k_2_subset_4000", 10, 3, 2),
        ("compas_3d_10k_2_subset_5000", 10, 3, 2),
        ("compas_3d_10k_2_subset_6889", 10, 3, 2),
        ("compas_3d_10k_2_subset_6889", 20, 3, 2),
        ("compas_3d_10k_2_subset_6889", 30, 3, 2),
        ("compas_3d_10k_2_subset_6889", 40, 3, 2),
        ("compas_3d_10k_2_subset_6889", 50, 3, 2),
    ]

    for exp_compas in exps_compas:
        protectedCompasRace, nonProtectedCompasRace, data, sort_time = compasData.createRace(exp_compas[0], exp_compas[2])
        rankAndDump(protectedCompasRace, nonProtectedCompasRace, exp_compas, pairsOfPAndAlpha, data, sort_time)


def createAndRankDOTData():

    pairsOfPAndAlpha = [(0.1, 0.0122),
                        (0.2, 0.0101),
                        (0.3, 0.0092),
                        (0.4, 0.0088),
                        (0.5, 0.0084),
                        (0.6, 0.0085),
                        (0.7, 0.0084),
                        (0.8, 0.0084),
                        (0.9, 0.0096)]

    exps_dot = [
        ("DOT_3d_1300k_rect_12_4_subset_10000", 10, 3, 12),
        ("DOT_3d_1300k_rect_12_4_subset_50000", 10, 3, 12),
        ("DOT_3d_1300k_rect_12_4_subset_100000", 10, 3, 12),
        ("DOT_3d_1300k_rect_12_4_subset_500000", 10, 3, 12),
        ("DOT_3d_1300k_rect_12_4_subset_1000000", 10, 3, 12),
    ]

    for exp_dot in exps_dot:
        protectedDotCarrier, nonProtectedDotCarrier, data, sort_time = dotData.createCarrier(exp_dot[0], exp_dot[2])
        rankAndDump(protectedDotCarrier, nonProtectedDotCarrier, exp_dot, pairsOfPAndAlpha, data, sort_time)


def createAndRankAllstarData():
    pairsOfPAndAlpha = [(0.1, 0.1),  # no real results, skip in evaluation
                        (0.2, 0.1),  # no real results, skip in evaluation
                        (0.3, 0.1),  # no real results, skip in evaluation
                        (0.4, 0.1),  # no real results, skip in evaluation
                        (0.5, 0.0168),
                        (0.6, 0.0321),
                        (0.7, 0.0293),
                        (0.8, 0.0328),
                        (0.9, 0.0375)]

    exps_allstar = [
        ("allstar_4d_1610_4_subset_400", 10, 4, 2),
        ("allstar_4d_1610_4_subset_800", 10, 4, 2),
        ("allstar_4d_1610_4_subset_1200", 10, 4, 2),
        ("allstar_4d_1610_4", 10, 4, 2),
    ]

    for exp_allstar in exps_allstar:
        protectedAllstarConf, nonProtectedAllstarConf, data, sort_time = allstarData.createConf(exp_allstar[0], exp_allstar[2])
        rankAndDump(protectedAllstarConf, nonProtectedAllstarConf, exp_allstar, pairsOfPAndAlpha, data, sort_time)


def createAndRankSyntheticData():

    pairsOfPAndAlpha = [(0.1, 0.0122),
                        (0.2, 0.0101),
                        (0.3, 0.0092),
                        (0.4, 0.0088),
                        (0.5, 0.0084),
                        (0.6, 0.0085),
                        (0.7, 0.0084),
                        (0.8, 0.0084),
                        (0.9, 0.0096)]

    exps_synthetic = [
        ("dataset_independent_10k_3d_rect_2", 10, 3, 2),
        ("dataset_independent_50k_3d_rect_2", 10, 3, 2),
        ("dataset_independent_100k_3d_rect_2", 10, 3, 2),
        ("dataset_independent_500k_3d_rect_2", 10, 3, 2),
        ("dataset_independent_1000k_3d_rect_2", 10, 3, 2),
        ("dataset_independent_100k_2d_rect_2", 10, 2, 2),
        ("dataset_independent_100k_4d_rect_2", 10, 4, 2),
        ("dataset_independent_100k_5d_rect_2", 10, 5, 2),
        ("dataset_anti_correlated_100k_3d_rect_2", 10, 3, 2),
        ("dataset_anti_correlated_100k_3d_rect_3", 10, 3, 3),
        ("dataset_anti_correlated_100k_3d_rect_4", 10, 3, 4),
        ("dataset_anti_correlated_100k_3d_rect_5", 10, 3, 5),
    ]

    for exp_synthetic in exps_synthetic:
        protectedSynthetic, nonProtectedSynthetic, data, sort_time = syntheticData.createSynthetic(exp_synthetic[0], exp_synthetic[2])
        rankAndDump(protectedSynthetic, nonProtectedSynthetic, exp_synthetic, pairsOfPAndAlpha, data, sort_time)


'''
def createAndRankGermanCreditData(k):

    pairsOfPAndAlpha = [(0.1, 0.1),  # no real results, skip in evaluation
                        (0.2, 0.1),  # no real results, skip in evaluation
                        (0.3, 0.0220),
                        (0.4, 0.0222),
                        (0.5, 0.0207),
                        (0.6, 0.0209),
                        (0.7, 0.0216),
                        (0.8, 0.0216),
                        (0.9, 0.0256)]

    protectedGermanCreditGender, nonProtectedGermanCreditGender = germanCreditData.create(
        "../rawData/GermanCredit/GermanCredit_sex.csv", "DurationMonth", "CreditAmount",
        "score", "sex", protectedAttribute=["female"])
    rankAndDump(protectedGermanCreditGender, nonProtectedGermanCreditGender, k,
                       "GermanCreditGender", "../results/rankingDumps/German Credit/Gender",
                       pairsOfPAndAlpha)

    protectedGermanCreditAge25, nonProtectedGermanCreditAge25 = germanCreditData.create(
        "../rawData/GermanCredit/GermanCredit_age25.csv", "DurationMonth", "CreditAmount",
        "score", "age25", protectedAttribute=["younger25"])
    rankAndDump(protectedGermanCreditAge25, nonProtectedGermanCreditAge25, k,
                       "GermanCreditAge25", "../results/rankingDumps/German Credit/Age25",
                       pairsOfPAndAlpha)

    protectedGermanCreditAge35, nonProtectedGermanCreditAge35 = germanCreditData.create(
        "../rawData/GermanCredit/GermanCredit_age35.csv", "DurationMonth", "CreditAmount",
        "score", "age35", protectedAttribute=["younger35"])
    rankAndDump(protectedGermanCreditAge35, nonProtectedGermanCreditAge35, k,
                       "GermanCreditAge35", "../results/rankingDumps/German Credit/Age35",
                       pairsOfPAndAlpha)


def createAndRankChileData(k):

    # loop through all files
    chileDir = '../rawData/ChileSAT/Dataset'
    pairsOfPAndAlpha = [(0.1, 0.0122),
                        (0.2, 0.0101),
                        (0.3, 0.0092),
                        (0.4, 0.0088),
                        (0.5, 0.0084),
                        (0.6, 0.0085),
                        (0.7, 0.0084),
                        (0.8, 0.0084),
                        (0.9, 0.0096)]
    for root, dirs, filenames in os.walk(chileDir):
        for chileFile in filenames:
            if not chileFile.startswith('.') and os.path.isfile(os.path.join(root, chileFile)):

                chileFile = '../rawData/ChileSAT/Dataset/' + chileFile
                print("reading: " + chileFile)

                protectedChileSATSchool, nonProtectedChileSATSchool = chileSATData.createSchool(chileFile, 6)
                rankAndDump(protectedChileSATSchool, nonProtectedChileSATSchool, k, "ChileSATSchool",
                                   "../results/rankingDumps/ChileSAT", pairsOfPAndAlpha)
                # protectedChileSATNat, nonProtectedChileSATNat = chileSATData.createNationality(chileFile, 6)
                # rankAndDump(protectedChileSATNat, nonProtectedChileSATNat, k, "ChileSATNationality",
                #                    "../results/rankingDumps/ChileSAT", pairsOfPAndAlpha)
'''


def i_dom_by_j(data_i, data_j, dim):
    le_count = 0
    l_count = 0
    for c in range(dim):
        if data_i[c] <= data_j[c]:
            le_count += 1
            if data_i[c] < data_j[c]:
                l_count += 1
    if le_count == dim and l_count > 0:
        return True
    else:
        return False


def getVioCount(list_data, fairRanking, fairNotSelected, k, dim):
    final_vio_count = 0

    for i in fairRanking:
        top_k_i = i.id
        vio_i = 0

        for j in fairNotSelected:
            remaining_j = j.id

            if i_dom_by_j(list_data[top_k_i - 1], list_data[remaining_j - 1], dim):
                vio_i += 1 ## violation if a top-k tuple is *dominated by* a remaining tuple

        final_vio_count += vio_i
        # print("vio_count for", i, ":", vio_i)

    vio = 0
    for i in range(k):
        data_i = list_data[fairRanking[i].id - 1]

        for j in range(i, k):
            data_j = list_data[fairRanking[j].id - 1]

            if i_dom_by_j(data_i, data_j, dim):
                vio += 1 ## violation if a toper tuple is *dominated by* a later tuple

    final_vio_count += vio
    # print("vio_count for top-k:", vio)

    return final_vio_count


def rankAndDump(protected, nonProtected, config, pairsOfPAndAlpha, data, sort_time=0):
    """
    creates all rankings we need for one experimental data set and writes them to disk to be used later

    @param protected:        list of protected candidates, assumed to satisfy in-group monotonicty
    @param nonProtected:     list of non-protected candidates, assumed to satisfy in-group monotonicty
    @param k:                length of the rankings we want to create
    @param dataSetName:      determines which data set is used in this experiment
    @param directory:        directory in which to store the rankings
    @param pairsOfPAndAlpha: contains the mapping of a certain alpha correction to be used for a certain p

    The experimental setting is as follows: for a given data set of protected and non-
    protected candidates we create the following rankings:
    * a colorblind ranking,
    * a ranking as in Feldman et al
    * ten rankings using our FairRankingCreator, with p varying from 0.1, 0.2 to 0.9, whereas alpha
      stays 0.1

    """

    dataSetName = config[0]
    k = config[1]
    dim = config[2]
    numGrp = config[3]

    directory = "-".join([dataSetName, str(k), str(dim), str(numGrp)])

    print("====================================================================")
    print("create rankings of {0}".format(dataSetName))

    result_folder = os.path.join(os.getcwd(), "..", "results", "results_fa-ir", directory)
    if not os.path.exists(result_folder):
        print("makedir:", result_folder)
        os.makedirs(result_folder)

    t1_start = process_time()

    print("fair rankings", end='', flush=True)

    run_count = 0
    fairSelecteds = []
    fairNotSelecteds = []
    for pair in pairsOfPAndAlpha:
        fairSelected, fairNotSelected = fairRanking(k, protected, nonProtected, pair[0], pair[1])
        fairSelecteds.append(fairSelected)
        fairNotSelecteds.append(fairNotSelected)
        run_count += 1

    print(" [Done]")

    t1_stop = process_time()

    list_data = data.values.tolist()

    # pdb.set_trace()

    total_time = sort_time + (t1_stop-t1_start)/10

    print("Sorting time:", sort_time)
    print("Elapsed time during the whole program in seconds:", total_time)

    for i in range(run_count):
        vio_count = getVioCount(list_data, fairSelecteds[i], fairNotSelecteds[i], k, dim)
        print("vio_count:", vio_count)

        for top_k_i in range(len(fairSelecteds[i])):
            top_k_id = fairSelecteds[i][top_k_i].id
            data_item = list_data[top_k_id - 1]
            print(top_k_i, ":", top_k_id, ",", data_item)

        result_txt = os.path.join(result_folder, "stat_fa_ir_" + str(i + 1) + ".txt")
        with open(result_txt, "w") as f:
            fairSelected = fairSelecteds[i]
            for top_k_i in range(len(fairSelected)):
                top_k_id = fairSelected[top_k_i].id
                f.write("{tkid} {tkgrp}\n".format(
                    tkid=top_k_id,
                    tkgrp=int(list_data[top_k_id - 1][len(list_data[top_k_id - 1]) - 1])
                ))
            f.write("{}\n".format(vio_count))

        time_txt = os.path.join(result_folder, "time_fa_ir_" + str(i + 1) + ".txt")
        with open(time_txt, "w") as f:
            f.write(str(total_time))


if __name__ == '__main__':
    main()
