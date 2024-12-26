# system lib
import argparse
import time
import warnings
import os

# tools
import torch
from torch.autograd import Variable
import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd
from functools import reduce
import magic
from tqdm import tqdm

# cscFEA lib
from ClassFlux import FLUX # Flux class network
from util import pearsonr
from DatasetFlux import MyDataset


# hyper parameters
LEARN_RATE = 0.008  
#EPOCH = 100
LAMB_BA = 1
LAMB_NG = 1
LAMB_CELL =  1
LAMB_RXN = 2
cellTypeSpecific = [[0.2, 2], [0.1, 4]]


# The and rule reflects that some reactions need subunits
# Or rule is the sum of gene expression levels, while and rule is the mean of those.
def andRule(mat, rules):
    mat['N/A'] = 0
    for i in range(mat.shape[0]):
        tmp = []
        tmptmp = []
        for _ in rules:
            try:
                tmptmp += _
            except:
                tmp.append('N/A')

        if tmptmp != ['']:
            tmptmp = list(set(tmptmp))
            andMean = ((mat.iloc[[i], :].loc[:, tmptmp].mean()).iloc[0])/len(tmptmp)
            mat.iloc[[i], :].loc[:, tmptmp] = andMean
            mat.loc[mat.index[i], [g for g in mat.columns if g not in tmptmp]] = 0
        else:
            mat.loc[mat.index[i], [g for g in mat.columns if g not in tmp]] = 0

    del mat['N/A']
    return mat


def myLoss(r, c, lamb1 = 0.2, lamb2= 0.2, lamb3 = 0.2, lamb4 = 0.2, geneScale = None, reactionScale = None, dirConstraints = None, balConstraints = None):    
    # balance constrain
    total1 = torch.pow(c, 2)
    total1 = torch.sum(total1, dim = 1) * balConstraints
    
    # non-negative constrain
    # Forward only reactions yield positive when estimated to be reverse, and vice versa.
    # Bidirectional reactions do not contribute to loss.
    error = torch.abs(r*dirConstraints) - (r*dirConstraints)
    total2 = torch.sum(error, dim=1) 
    
    # sample-wise variation constrain 
    diff = torch.pow(torch.sum(r, dim=1) - geneScale, 2)
    if sum(diff > 0) == r.shape[0]: # solve Nan after several iteraions
        total3 = torch.pow(diff, 0.5)
    else:
        total3 = diff
    
    # reaction-wise variation constrain
    if lamb4 > 0 :
        corr = torch.FloatTensor(np.ones(r.shape[0]))
        for i in range(r.shape[0]):
            # if pearson r cannot be calculated
            # (when no gene associated with any reaction is expressed)
            # exclude from loss function.
            if reactionScale[i, :].sum() == 0:
                corr[i] = 1
            else:
                corr[i] = pearsonr(r[i], reactionScale[i])
        corr = torch.abs(corr)
        penal_r_var = torch.FloatTensor(np.ones(r.shape[0])) - corr
        total4 = penal_r_var
    else:
        total4 = torch.FloatTensor(np.zeros(r.shape[0]))
            
    # loss
    loss1 = torch.sum(lamb1 * total1)
    loss2 = torch.sum(lamb2 * total2)
    loss3 = torch.sum(lamb3 * total3)
    loss4 = torch.sum(lamb4 * total4)

    loss = loss1 + loss2 + loss3 + loss4
    return loss, loss1, loss2, loss3, loss4
    

def main(args):

    # set arguments
    data_path = args.data_dir
    input_path = args.input_dir
    res_dir = args.res_dir
    test_file = args.test_file
    cm_file = args.stoichiometry_matrix
    sc_imputation = args.sc_imputation
    cName_file = args.cName_file
    fileName = args.output_flux_file
    balanceName = args.output_balance_file
    EPOCH = args.train_epoch
    gene_reaction_rule_path = args.gene_reaction_rule_path
    reaction_direction_path = args.reaction_direction_path
    scoring = args.scoring

    if EPOCH <= 0:
        raise NameError('EPOCH must greater than 1!')

    # choose cpu or gpu automatically
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    # read data
    print("Starting load data...")

    # read table of gene expression level
    geneExpr = pd.read_csv(
                input_path + '/' + test_file,
                index_col=0)

    # cellType Specificity :
    # for the top 30% of genes of each cell or cluster, multiply 4 with gene expression level
    # for the top 30-40% of genes of them, multiply 3 with gene expression level for default

    print('Applying cellType specificity')
    weightAlign = np.ones(geneExpr.shape[0])

    for _ in cellTypeSpecific:
        weightAlign[int(-(_[0]*len(weightAlign))): ] = _[1]

    for topCombi in range(geneExpr.shape[1]):
        singleColumn = np.array(geneExpr[geneExpr.columns[topCombi]].values)
        singleColRank = singleColumn.argsort()
        singleColWeight = np.ones(len(singleColumn))

        for alignCount in range(len(weightAlign)):
            singleColWeight[singleColRank[alignCount]] = weightAlign[alignCount]

        if topCombi == 0:
            geneExprVal = pd.DataFrame(singleColumn*singleColWeight)
        else:
            geneExprVal = pd.concat([geneExprVal, pd.DataFrame(singleColumn*singleColWeight)], axis = 1)

    geneExprVal.index = geneExpr.index
    geneExprVal.columns = geneExpr.columns
    geneExpr.iloc[:, :] = geneExprVal.iloc[:, :]
    del geneExprVal

    geneExpr = geneExpr.T
    geneExpr = geneExpr * 1.0

    if sc_imputation == True:
        magic_operator = magic.MAGIC()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            geneExpr = magic_operator.fit_transform(geneExpr)

    if geneExpr.max().max() > 50:
        geneExpr = (geneExpr + 1).apply(np.log2)  

    geneExprSum = geneExpr.sum(axis=1)
    stand = geneExprSum.mean()
    geneExprScale = geneExprSum / stand
    geneExprScale = torch.FloatTensor(geneExprScale.values).to(device)
    
    BATCH_SIZE = geneExpr.shape[0]

    
    # load gene reaction rule
    with open(data_path + '/' + gene_reaction_rule_path, 'r') as grrF:
        grr = grrF.read()
    grr = grr.split('\n')

    if grr[-1] == '':
        del grr[-1]
    else:
        pass

    # sort gene reaction rules by reaction name
    rxnDict = {}
    for  _ in grr:
        try:
            rxnDict[_.split(',')[0]].append(_.split(',')[1].split(' and '))
        except:
            rxnDict[_.split(',')[0]] = [_.split(',')[1].split(' and ')]

    # set of all gene included in reactions
    tmpreactionGene = [_.split(',')[1].split(' and ') for _ in grr]
    reaction_gene_all = set(reduce(lambda grrTmpA, grrTmpB: grrTmpA + grrTmpB, tmpreactionGene, []))
    del tmpreactionGene


    # fix genes in both data and reaction rules
    # set virtual genes in reaction rules but not in data
    data_gene_all = set(geneExpr.columns)
    gene_overlap = list(data_gene_all.intersection(reaction_gene_all))
    virtual_gene = list(set(reaction_gene_all) - set(data_gene_all) - set(['']))
    gene_overlap.sort()

    cmMat = pd.read_csv(
            data_path + '/' + cm_file,
            sep=',',
            header=None)
    cmMat = cmMat.values
    cmMat = torch.FloatTensor(cmMat).to(device)
    
    if cName_file != 'noCompoundName':
        print("Load compound name file, the balance output will have compound name.")
        cName = pd.read_csv(
                data_path + '/' + cName_file,
                sep=',',
                header=0)
        cName = cName.columns
    print("Load data done.")
    
    print("Starting process data...")

    # extract overlap gene
    geneExpr = geneExpr[gene_overlap]
    gene_names = geneExpr.columns
    cell_names = geneExpr.index.astype(str)
    geneExprIndex = geneExpr.index.tolist()

    n_rxn = len(rxnDict)
    n_genes = len(gene_names)
    n_cells = len(cell_names)
    n_comps = cmMat.shape[0]
    geneExprDf = pd.DataFrame(columns = ['reaction_gene'] + list(cell_names))

    
    rxnRules = list(rxnDict.values())
    for i in range(n_rxn):
        temp = andRule(geneExpr.copy(), rxnRules[i])
        temp = temp.T
        temp['reaction_gene'] = ['%02d_%s' % (i,g) for g in gene_names]
        geneExprDf = pd.concat([geneExprDf, temp], ignore_index = True)


    # Set Direction Constraints.
    # Reactions are grouped into forward only, bidirectional and reverse only.
    # Groups are based on Recon3D model.
    reaction_direction = pd.read_csv(data_path + '/' + reaction_direction_path, index_col = 0)
    rdDict = {}
    for i in range(reaction_direction.shape[0]):
        rdDict[reaction_direction.iloc[i, 0]] = reaction_direction.iloc[i, 1]
    
    dirCons = []
    for _ in rxnDict.keys():
        dirCons.append(rdDict[_])

    dirCons = torch.FloatTensor(np.array(dirCons))
    dirFactor = (dirCons+0.1)/abs(dirCons+0.1)
    dirCons = torch.pow(dirCons, 2)
    cmMat = cmMat*dirFactor

    geneExprDf.index = geneExprDf['reaction_gene']
    geneExprDf.drop('reaction_gene', axis = 'columns', inplace = True)
    X = geneExprDf.values.T.astype(np.float32)
    X = torch.FloatTensor(X).to(device)

    balanceCons = torch.where(torch.where(cmMat != 0, 1, 0).sum(axis=1) > 1, 1, 0)
    balanceCons = balanceCons.reshape(len(balanceCons), 1)
    
    # prepare data for constraint of reaction variation based on gene
    df = geneExprDf
    df.index = [i.split('_')[0] for i in df.index]
    df.index = df.index.astype(int)   # mush change type to ensure correct order, T column name order change!
    rxn_scale = df.groupby(df.index).sum().T
    rxn_scale = torch.FloatTensor(rxn_scale.values)

    print("Process data done.")
    
# =============================================================================
    #NN
    torch.manual_seed(3686)
    net = FLUX(X, n_rxn, f_in = n_genes, f_out = 1).to(device)
    optimizer = torch.optim.Adam(net.parameters(), lr = LEARN_RATE)

    #Dataloader
    dataloader_params = {'batch_size': BATCH_SIZE,
                         'shuffle': False,
                         'num_workers': 0,
                         'pin_memory': False}

    dataSet = MyDataset(X, geneExprScale, rxn_scale)
    train_loader = torch.utils.data.DataLoader(dataset=dataSet,
                                               **dataloader_params)
    
# =============================================================================

  
   
    
# =============================================================================
    print("Starting train neural network...")
    start = time.time()  
#   training
    loss_v = []
    loss_v1 = []
    loss_v2 = []
    loss_v3 = []
    loss_v4 = []
    net.train()
    timestr = time.strftime("%Y%m%d-%H%M%S")
    lossName = "./output/lossValue_" + timestr + ".txt"
    file_loss = open(lossName, "a")
    for epoch in tqdm(range(EPOCH)):
        loss, loss1, loss2, loss3, loss4 = 0,0,0,0,0
        
        for i, (X, X_scale, r_scale) in enumerate(train_loader):

            X_batch = Variable(X.float().to(device))
            X_scale_batch = Variable(X_scale.float().to(device))
            r_scale_batch = Variable(r_scale.float().to(device))
            
            out_r_batch, out_c_batch = net(X_batch, n_rxn, n_genes, n_comps, cmMat)
            loss_batch, loss1_batch, loss2_batch, loss3_batch, loss4_batch  = myLoss(out_r_batch, out_c_batch, 
                                                                                         lamb1 = LAMB_BA, lamb2 = LAMB_NG, lamb3 = LAMB_CELL, lamb4 = LAMB_RXN, 
                                                                                         geneScale = X_scale_batch, reactionScale = r_scale_batch, dirConstraints = dirCons, balConstraints = balanceCons)
            
            optimizer.zero_grad()
            loss_batch.backward()
            optimizer.step()
            
            loss += loss_batch.cpu().data.numpy()
            loss1 += loss1_batch.cpu().data.numpy()
            loss2 += loss2_batch.cpu().data.numpy()
            loss3 += loss3_batch.cpu().data.numpy()
            loss4 += loss4_batch.cpu().data.numpy()
            
        #print('epoch: %02d, loss1: %.8f, loss2: %.8f, loss3: %.8f, loss4: %.8f, loss: %.8f' % (epoch+1, loss1, loss2, loss3, loss4, loss))
        file_loss.write('epoch: %02d, loss1: %.8f, loss2: %.8f, loss3: %.8f, loss4: %.8f, loss: %.8f. \n' % (epoch+1, loss1, loss2, loss3, loss4, loss))
        
        loss_v.append(loss)
        loss_v1.append(loss1)
        loss_v2.append(loss2)
        loss_v3.append(loss3)
        loss_v4.append(loss4)
        #pd.concat(loss_v, loss, axis = 0)
        #pd.concat(loss_v1, loss1, axis = 0)
        #pd.concat(loss_v2, loss2, axis = 0)
        #pd.concat(loss_v3, loss3, axis = 0)
        #pd.concat(loss_v4, loss4, axis = 0)
# =============================================================================
    end = time.time()
    print("Training time: ", end - start) 
    
    file_loss.close()
    plt.plot(loss_v, '--')
    plt.plot(loss_v1)
    plt.plot(loss_v2)
    plt.plot(loss_v3)
    plt.plot(loss_v4)
    plt.legend(['total', 'balance', 'negative', 'cellVar', 'reactionVar']);
    imgName = './' + res_dir + '/loss_' + timestr + ".png"
    plt.savefig(imgName)
    plt.clf()
    timeName =  './' + res_dir + '/time_' + timestr + ".txt"
    f = open(timeName, "a")
    runTimeStr = str(end - start)
    f.write(runTimeStr)
    f.close()   

    
#    Dataloader
    dataloader_params = {'batch_size': 1,
                         'shuffle': False,
                         'num_workers': 0,
                         'pin_memory': False}

    dataSet = MyDataset(X, geneExprScale, rxn_scale)
    test_loader = torch.utils.data.DataLoader(dataset=dataSet,
                          **dataloader_params)
   
    #testing
    fluxStatuTest = np.zeros((n_cells, n_rxn), dtype='f') #float32
    balanceStatus = np.zeros((n_cells, n_comps), dtype='f')
    net.eval()
    for epoch in range(1):
        loss, loss1, loss2 = 0,0,0
        
        for i, (X, X_scale, _) in enumerate(test_loader):

            X_batch = Variable(X.float().to(device))
            out_r_batch, out_c_batch = net(X_batch, n_rxn, n_genes, n_comps, cmMat)
            
            # save data
            fluxStatuTest[i, :] = out_r_batch.detach().numpy()
            balanceStatus[i, :] = out_c_batch.detach().numpy()
            
                   
    
    # save to file
    if fileName == 'NULL':
        # user do not define file name of flux
        if '/' in test_file:
            fileName = "./" + res_dir + "/" + test_file.split('/')[-1][0:-4] + "_reaction" + str(n_rxn) + "_cell" + str(n_cells) + "_batch" + str(BATCH_SIZE) + \
                        "_LR" + str(LEARN_RATE) + "_epoch" + str(EPOCH) + "_SCimpute_" + str(sc_imputation)[0] + \
                        "_lambBal" + str(LAMB_BA) + "_lambSca" + str(LAMB_NG) + "_lambCellCor" + str(LAMB_CELL) + "_lambRxn" + str(LAMB_RXN) + \
                        '_' + timestr + ".csv"
        else:
            fileName = "./" + res_dir + "/" + test_file[-len(test_file):-4] + "_reaction" + str(n_rxn) + "_cell" + str(n_cells) + "_batch" + str(BATCH_SIZE) + \
                        "_LR" + str(LEARN_RATE) + "_epoch" + str(EPOCH) + "_SCimpute_" + str(sc_imputation)[0] + \
                        "_lambBal" + str(LAMB_BA) + "_lambSca" + str(LAMB_NG) + "_lambCellCor" + str(LAMB_CELL) + "_lambRxn" + str(LAMB_RXN) + \
                        '_' + timestr + ".csv"

    setF = pd.DataFrame(fluxStatuTest)
    try:
        setF.columns = list(rxnDict.keys())
    except:
        setF.columns = geneList

    setF.index = geneExprIndex
    setF.to_csv(fileName)

    
    setB = pd.DataFrame(balanceStatus)
    setB.rename(columns = lambda x: x + 1)
    setB.index = setF.index
    if cName_file != 'noCompoundName':
        setB.columns = cName
    if balanceName == 'NULL':
        # user do not define file name of balance
        balanceName = "./output/balance_" + timestr + ".csv"
    setB.to_csv(balanceName)


    # for visualization, start further process
    # it includes calculation of mean values by groups, min-max normalization, and scoring.
    # pairs of metabolites where the disapperance of one corresponds to the formation of the other in each reaction are picked.
    # based on the sequence of these conversions, a hierarchical structure of each metabolic reaction is established.

    if scoring:
        # set groups for comparison, and set hierarchical structure of reactions
        compare_groups = [['1', '2', '3'], ['4'], ['5'], ['6']]
        group_names = ['Preflight', 'R+1', 'R+45', 'R+82']

        if len(compare_groups) != len(group_names):
            print('Number of names for comparison groups are not matched. Rename as compare_groups')
        reaction_hierarchy = [['NADH2_u10mi', 'FADH2ETC', 'r0205'], ['CYOR_u10mi'], ['CYOOm2i', 'CYOOm3i']]
        reformat_table = setF.copy()
        group_sep = []
        for _ in list(reformat_table.index):
            sep_tmp = str(_.split('_')[1])
            for groups in compare_groups:
                if sep_tmp in groups:
                    group_sep.append(groups[0])

        # calculate mean flux level
        reformat_table['groups'] = group_sep
        reformat_table_mean_values = reformat_table.groupby('groups').mean()
        reformat_table_mean_values = np.array(reformat_table_mean_values)

        # min-max normalization
        val_max = reformat_table_mean_values.max().max()
        val_min = reformat_table_mean_values.min().min()

        normalized_mv = (reformat_table_mean_values - val_min) / (val_max - val_min)

        scores = np.zeros(normalized_mv.shape)
        mins = np.ones(normalized_mv.shape)
        meanRef = normalized_mv.mean(axis = 0)
        reactionOrderDict = {list(setF.columns)[i] : i for i in range(len(setF.columns))}
        minThreshold = 0.05
        leftReactions = list(setF.columns)

        # Scoring
        # reactions coming earlier in the conversion of pairs were assigned higher score.
        # for simultaneous reactions, higher scores were assigned based on relative mean flux level (MFL).
        # among same reactions, samples with the lowest MFL among same reactions were assigned lower score.
        # reactions with too low MFL were assigned 0, except for the highest MFL in each reaction.
        # the other reactions that do not involve in the conversion of pairs were assigned the same score of the higest assumed to be maximized.
        
        for reaction_order in reaction_hierarchy:
            asNum = []
            for _ in reaction_order:
                asNum.append(reactionOrderDict[_])
                leftReactions.remove(_)

            if len(asNum) >= 3:
                tmpMax = asNum[meanRef[asNum].argmax()]
                tmpMin = asNum[meanRef[asNum].argmin()]

                scores[:, tmpMax] += 2
                scores[normalized_mv[:, tmpMax].argmin(), tmpMax] += -2
                
                if meanRef[tmpMin] > minThreshold:
                    scores[:, tmpMin] += -2
                else:
                    mins[normalized_mv[:, tmpMin].argmax(), tmpMin] += 1
                    mins[:, tmpMin] += -1

                for i in asNum:
                    if i == tmpMax:
                        pass
                    elif i == tmpMin:
                        pass
                    else:
                        tmp = np.argsort(normalized_mv[:, i])
                        scores[tmp[int(len(tmp)/2):], i] += 2
                        scores[tmp[:int(len(tmp)/2)], i] += 1

            elif len(asNum) == 2:
                for i in range(2):
                    if meanRef[asNum[i]] > minThreshold:
                        scores[normalized_mv[:, asNum[i]].argmin(), asNum[i]] += -2
                    else:
                        mins[normalized_mv[:, asNum[i]].argmax(), asNum[i]] += 1
                        mins[:, asNum[i]] += -1

            else:
                if meanRef[asNum[0]] > minThreshold:
                    scores[normalized_mv[:, asNum].argmin(), asNum[0]] += -2
                else:
                    mins[normalized_mv[:, asNum[0]].argmax(), asNum[0]] += 1
                    mins[:, asNum[0]] += -1

            scores += 1
        
        tmp = []
        for _ in leftReactions:
            tmp.append(reactionOrderDict[_])
        scores[:, tmp] = scores.max().max()

        scores = (abs(scores) + scores)/2
        scores = scores * mins

        score_df = pd.DataFrame(scores)
        score_df.columns = setF.columns
        score_df.index = group_names
        score_df.to_csv(f'{fileName[:-4]}_score.csv')

        

    print("cscFEA job finished. Check result in the desired output folder.")
    
    
    return


def parse_arguments(parser):


    parser.add_argument('--data_dir', type=str, default='data', metavar='<data_directory>',
                        help='The data directory for cscFEA model files.')
    parser.add_argument('--input_dir', type=str, default='input', metavar='<input_directory>',
                        help='The data directory for single cell input data.')
    parser.add_argument('--res_dir', type=str, default='output', metavar='<data_directory>',
                        help='The data directory for result [output]. The output of cscFEA includes two matrices, predicted metabolic flux and metabolites stress at single cell or single cluster resolution.')
    parser.add_argument('--test_file', type=str, default='insp4_mean_merged.csv', 
                        help='The test SC file [input]. The input of cscFEA is a single cell profile matrix, where row is gene and column is cell. An example dataset is provided in /data/ folder. The input can be raw counts or normalised counts. The logarithm would be performed if value larger than 50.')
    parser.add_argument('--stoichiometry_matrix', type=str, default='cmMat_Recon3D_OXPHOS.csv', 
                        help='The table describes relationship between compounds and reactions. Each row is an intermediate metabolite and each column is metabolic reaction. Stoichiometry matrix of OXPHOS from Recon3D is provided in /data/ folder.')
    parser.add_argument('--cName_file', type=str, default='cName_Recon3D_OXPHOS.csv',
                        help='The name of compounds. The table contains two rows. First row is compounds name and second row is corresponding id.')
    parser.add_argument('--sc_imputation', type=eval, default='False', choices=[True, False],
                        help='Whether perform imputation for SC dataset.')
    parser.add_argument('--output_flux_file', type=str, default='NULL', 
                        help='User defined predicted flux file name.')
    parser.add_argument('--output_balance_file', type=str, default='NULL', 
                        help='User defined predicted balance file name.')
    parser.add_argument('--train_epoch', type=int, default=100, nargs='?',
                        help='User defined EPOCH (training iteration).')
    parser.add_argument('--gene_reaction_rule_path', type=str, default='geneReactionRule_Recon3D_OXPHOS.csv', metavar='<data_directory>',
                        help='The table of gene reaction rules. It contains two columns. The first column is name of reaction, and the second column is corresponding reaction formula. Genes linked with OR are separated into different rows with same name. Genes linke with AND are written, e.g. A and B and C. For the details, please check enclosed sample.')
    parser.add_argument('--reaction_direction_path', type=str, default='reaction_direction_Recon3D.csv', metavar='<data_directory>',
                        help='The table of directionalities of reactions. It contains two columns. The first column is name of reaction, and the second column is corresponding directionality based on Recon3D. 1 for forward only, 0 for bidirectional, -1 for backward reaction only.')
    parser.add_argument('--scoring', type=eval, default=False, choices=[True, False],
                        help='To ease visualizing relative differences between groups, grouping, calculating mean, normalizing, and scoring are introduced.')
    
    
    

    args = parser.parse_args()

    return args


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='cscFEA: A graph neural network model to estimate cell type-specific cell-wise metabolic flux using single cell RNA-seq data')
    args = parse_arguments(parser)
    main(args)

