
# Libraries
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import time


    # Fourth step
# Volcano plot



# This file contain the names of the clusters (columsn on the snRNA-seq) sort as snRNA-seq matrix
iterar=pd.read_csv("fenotipos.csv", index_col=0)
iterar=iterar["0"].tolist()
print("\n Fenotipos \n")
print(iterar)



for phenotype in iterar:
   
    gene_expr=pd.read_csv(phenotype+"_pvalue_foldchange.csv", index_col=0)
    #gene_expr=gene_expr[["p_values","fold_change"]]
    print(gene_expr)


    plt.scatter(x=gene_expr['fold_change'],y=gene_expr['p_values'].apply(lambda x:-np.log10(x)),s=1,label="Not significant")

    # highlight down- or up- regulated genes
    down = gene_expr[(gene_expr['fold_change']<=-2)&(gene_expr['p_values']<=0.001)] #Only genes with p_value less than or equal to 0.001 will be plotted for down genes.
    up = gene_expr[(gene_expr['fold_change']>=2)&(gene_expr['p_values']<=0.001)] #Only genes with p_value less than or equal to 0.001 will be plotted for up genes.

    #Plot the data
    plt.scatter(x=down['fold_change'],y=down['p_values'].apply(lambda x:-np.log10(x)),s=3,label="Down-regulated",color="blue")
    plt.scatter(x=up['fold_change'],y=up['p_values'].apply(lambda x:-np.log10(x)),s=3,label="Up-regulated",color="red")

    plt.xlabel("fold_change")
    plt.ylabel("p_value")
    plt.axvline(-2,color="grey",linestyle="--")
    plt.axvline(2,color="grey",linestyle="--")
    plt.axhline(2,color="grey",linestyle="--")
    plt.legend()


    print("Genes up-regulated \n")
    print(up)
    # Save the up regulated genes
    up.to_csv(phenotype+"_up.csv")


    print("Genes down-regulated \n")
    print(down)
    #Save the down regulared genes
    down.to_csv(phenotype+"_down.csv")

    # Save the volcano plot figure
    plt.savefig(phenotype+"_volcano_plot.jpg")
    plt.show()

print("Finished")

