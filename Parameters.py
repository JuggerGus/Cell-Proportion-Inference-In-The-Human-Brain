# Libraries
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import time




# Read the snRNA-seq sort by phenotype (Ej. GABA_1, GABA_1....Gluata_1, Gluta_2) or by Class.
# The rows are genes and the columns cells
sort=pd.read_csv("merged_matrix_transpose.csv", index_col=0, nrows=10)
print("\n Gene expression matrix \n")
print(sort)


# In case the snRNA-seq array has phenotype (clusters) names, a file with the range of phenotypes belonging to the main class (e.g. GABAergic, Glutamatergic, Non-neuronal) is needed.
range_types=pd.read_csv("cluster_lables_sort.csv")
print("\n Cluster ranges \n")
print(range_types)


# This file contain the names of the clusters (columsn on the snRNA-seq) sort as snRNA-seq matrix
iterar=pd.read_csv("fenotipos.csv", index_col=0)
iterar=iterar["0"].tolist()
print("\n Fenotipos \n")
print(iterar)




# p_value and Fold Change parameters for each gene in the each group (GABAergic, Glutamatergic, Non-neuronal)

for phenotype in iterar:

        # In this section we sort each cell belonging to main class (GABAergic, Glutamatergic, Non-neuronal) in the first iteration
        # the "cell" variable contain all cell types that belonging to class of interest to have the parameters (p_value and Fold), and the "noncell"
        # variable has the remaining cells. In this first iteration we obtain the parameters for each gene in the class of interest, and this cycle will repeat
        # this cycle will be repeated until all the groups are completed (the main cell classes).  


    # Obtain indexes to use as ranges
    cell_type=range_types["Cell type"]==phenotype
    print("\n Para cada phenotype obtenemos los TRUE \n" ) #Check point
    print(cell_type)


    # Create range for specific cell type "group1" change cell type of interest
    cell=range_types.loc[range_types["Cell type"]==phenotype]
    cell=cell["Unnamed: 0"].to_list()
    # Range of phenotype of interest
    print(cell)

    # Create range for the other cell types "group2" change cell type of interest
    noncell=range_types.loc[range_types["Cell type"]!=phenotype]
    noncell=noncell["Unnamed: 0"].to_list()
    # Range of the rest of the phenotypes
    print(noncell)



        # First step, separate the matrix into comparison groups.

    # Define two groups of samples (e.g. Neruon specific type vs. all others)
    # Depending on the ranks obtained above, the ranks to be taken for each of the groups will be as follows

    group1 = sort.columns[cell]  # Phenotype to compare with all others 0-xxx
    group2 = sort.columns[noncell] # Rest of phenotypes xxx-xxx-

    # Data visualization

    print("\n View the groups \n" ) #Check point
    print("Neruon specific type \n")
    print(group1)
    print("Others \n")
    print(group2)



        # Second step

    # Calculate the p-value for each gene using a t-test


    print("\n P_values \n" ) #Check point
    p_values = []
    #t_statitic= []
    # For each gene, the data for that gene is compared with all the columns
    for gene in sort.index:
        t, p = ttest_ind(sort.loc[gene, group1], sort.loc[gene, group2], equal_var=False)
        p_values.append(p)

        print("View p_value in each group",sort.loc[gene, group1])
        sort.loc[gene, group1].to_csv("Dataframe_group1.csv")
        sort.loc[gene, group2].to_csv("Dataframe_group2.csv")
        

        #t_statitic.append(t)

    print("P_values \n" ) #Check point
    print(p_values)


    sort["p_values"]=p_values
    #sort["t_statitic"]=t_statitic


        #Third step
    
    # Calculate the fold change for each gene
    print("Fold change \n" ) #Check point
    fold_changes = (sort[group2].mean(axis=1) / sort[group1].mean(axis=1)).apply(pd.np.log2)
    print("Fold change \n")
    print(fold_changes)

    print("\n Fold change_grupo2 \n")
    print(sort[group2].mean(axis=1))

    print("\n Fold change_grupo1 \n")
    print(sort[group1].mean(axis=1))


    # Add the fold changes to the gene expression matrix
    sort['fold_change'] = fold_changes
    #print(sort)


    # Save datafrmar with parameter for each gene of the class of interest 
    sort[["p_values","fold_change"]].to_csv(phenotype+"_pvalue_foldchange.csv")


    print("El proceso termino en: " +phenotype)
    #time.sleep(60)


print("Finished")
#time.sleep(60)

