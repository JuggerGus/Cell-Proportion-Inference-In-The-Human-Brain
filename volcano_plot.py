#LIbrerias
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import time


# Path para guardar el analisis por cluster
path="/home/r2d2/Desktop/Gus/Data_allen/" #Matrix
path_range="/home/r2d2/Desktop/Gus/Data_allen/Volcano/" #Clusters que sirven como header
path_volcano="/home/r2d2/Desktop/Gus/Data_allen/Volcano/Volcano_plot_per_phenotype/" #DataFrame con los p_values y fold_change
path_pvfc="/home/r2d2/Desktop/Gus/Data_allen/Volcano/Volcano_plot/" #Resultados, grafica volcano, down y up genes 


# Path de prueba, la carptea se llama "Carpeta_de_prueba"
path_prueba="/home/r2d2/Desktop/Gus/Data_allen/Volcano/Carpeta_de_prueba/"

# Path para guardar el analisis por clases 
path_class_pvfc="/home/r2d2/Desktop/Gus/Data_allen/Volcano/Analisis_clase_pvfc/"
path_class_down_up="/home/r2d2/Desktop/Gus/Data_allen/Volcano/Analisis_clase_up_down/"



# Archivo segun el analisis (clase o fenotipo)
    #Para Fenotipo: cluster_labels_sort.csv
    #Para Clase: Class_sort.csv
range_types=pd.read_csv(path_range+"Class_sort.csv")
print("\n Cluster ranges \n")
range_types.rename(columns={"0":"Cell type"}, inplace=True)
print(range_types)



# Lectura de archivo principal (matrix de expresion)
sort=pd.read_csv("sort_class.csv", index_col=0)
print("\n Gene expression matrix \n")
print(sort,"\n")
#time.sleep(100)

#new_colums=range_types["Cell type"].tolist()
#sort.columns=range_types["Cell type"]
#print(sort.columns)
#print(sort,"\n")


# Archivo segun el analisis (clase o fenotipo)
    #Para Fenotipo: fenotipos.csv
    #Para Clase: Class_fenotipo.csv
iterar=pd.read_csv("Class_fenotipo.csv", index_col=0)
iterar=iterar["0"].tolist()
print("\n Fenotipos \n")
print(iterar)




#Ciclo para obtener p-value y fold change para cada tipo celuar en cada vuelta del ciclo

for phenotype in iterar:
    #-----------------Primer paso-----------------
    # separar los grupos de comparacion (grupo 1 y grupo 2)


    # Obtener los indices para utilizarlos como rangos 
    cell_type=range_types["Cell type"]==phenotype
    print("\n Para cada phenotype obtenemos los TRUE \n" ) #Check point
    print(cell_type)


    # Creacion de rango para el tipo celualr especifico "group1" cambiar el tipo celular de interes
    cell=range_types.loc[range_types["Cell type"]==phenotype]
    cell=cell["Unnamed: 0"].to_list()
    # Rango del fenotipo ce interes
    print("\nRango del grupo 1 \n")
    print(cell)
    #time.sleep(30)


    # Creacion de rango para el los otrs tipos celulares "group2" cambiar el typo celular de interes 
    noncell=range_types.loc[range_types["Cell type"]!=phenotype]
    noncell=noncell["Unnamed: 0"].to_list()
    # Rango del resto de los fenotipos
    print("\nRango del grupo 2 \n")
    print(noncell)
    #time.sleep(30)


    # Define two groups of samples (e.g. Neruon specific type vs. all others)
    # Dependiendo de los rangos obtenidos arriba, seran los rangos a tomar para cada unos de los grupos
    group1 = sort.columns[cell]  # Fenotipo a comparar con todos los demas 0-xxx
    group2 = sort.columns[noncell] # Resto de los fenotipos xxx-xxx-

    # Vizualizacion de datos
    print("\nVer los dos grupos \n" ) #Check point
    print("Neruon specific type \n")
    print(group1)
    print("Others \n")
    print(group2)
    #time.sleep(100)
    
    # 20 minutes aprox


    #-----------------Segundo paso-----------------
    # Calculate the p-value for each gene using a t-test


    print("\n Entrada a proceso de p_values \n" ) #Check point
    p_values = []
    t_statitic= []
    # Para cada gen se comparan los datos de ese gen con todas las columnas 
    #grupos normales grupo 1 vs grupo 2
    for gene in sort.index:
        t, p = ttest_ind(sort.loc[gene, group2], sort.loc[gene, group1], equal_var=False)
        p_values.append(p)
        t_statitic.append(t)

        #print("\nHERE\n")
        print("\nGroup 1 \n",sort.loc[gene, group1])
        print("\nGroup 2 \n",sort.loc[gene, group2])
        #time.sleep(10)
        #sort.loc[gene, group1].to_csv(path_prueba+"Datafrmae_prueba_group1.csv")
        #sort.loc[gene, group2].to_csv(path_prueba+"Datafrmae_prueba_group2.csv")
        

        #t_statitic.append(t)

    print("P_values \n" ) #Check point
    print(p_values)

    sort["p_values"]=p_values
    sort["t_statitic"]=t_statitic
    #p_values=sort["p_values"]
    #print(sort)

    # Cambiar 
    #p_values.to_csv(path_volcano+"Exc L2-3 LINC00507 RPL9P17_pvalue.csv")


    #-----------------Tercer paso-----------------
    # Calculate the fold change for each gene


    print("Entrada a proceso de fold change \n" ) #Check point
    #grupos normales grupo 2 vs grupo 1
    fold_changes = (sort[group1].mean(axis=1) / sort[group2].mean(axis=1)).apply(pd.np.log2)
    print("Fold change \n")
    print(fold_changes)

    print("\n Media grupo 2 \n")
    print(sort[group2].mean(axis=1))

    print("\n Medioa grupo 1 \n")
    print(sort[group1].mean(axis=1))

    # Add the fold changes to the gene expression matrix
    sort['fold_change'] = fold_changes
    #print(sort)

    # CAMBIAR LABEL (Para la prueba se cambio el path, pero el original es "path_volcano")
        #Para Cluster: path_volcano
        #Para Clase: path_class_pvfc
        #Para prueba: path_prueba
    sort[["p_values","fold_change"]].to_csv(path_class_pvfc+phenotype+"_pvalue_foldchange_grupos_cambiados.csv")


    print("\nEl proceso va en: " +phenotype +"\n")
    print(sort[["p_values","t_statitic","fold_change"]])
    #time.sleep(1000)


print("Ya estuvo compa")
#time.sleep(10)





#25 minutos para terminar para ccada fenotipo aprox, para los 120 son aprox 3 dias




"""
#-----------------Cuarto paso-----------------
# Volcano plot


for phenotype in iterar:
   
    # CAMBIAR LABEL
    gene_expr=pd.read_csv(path_class_pvfc+phenotype+"_pvalue_foldchange_grupos_cambiados.csv", index_col=0)
    #gene_expr=gene_expr[["p_values","fold_change"]]
    print(gene_expr)


    plt.scatter(x=gene_expr['fold_change'],y=gene_expr['p_values'].apply(lambda x:-np.log10(x)),s=1,label="Not significant")

    # highlight down- or up- regulated genes
    down = gene_expr[(gene_expr['fold_change']<=-2)&(gene_expr['p_values']<=0.001)] #Solo los genes con p_value menor o igual a 0.001 se graficaran para genes down
    up = gene_expr[(gene_expr['fold_change']>=2)&(gene_expr['p_values']<=0.001)] #Solo los genes con p_value menor o igual a 0.001 se graficaran para genes up


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
    


    # CAMBIAR LABEL
        #Para Cluster: path_pvfc
        #Para Clase: path_class_down_up
    up.to_csv(path_class_down_up+phenotype+"_up_grupos_cambiados.csv")


    print("Genes down-regulated \n")
    print(down)

    # CAMBIAR LABEL
    down.to_csv(path_class_down_up+phenotype+"_down_grupos_cambiados.csv")

    # CAMBIAR LABEL
    plt.savefig(path_class_down_up+phenotype+"_volcano_plot_grupos_cambiados.jpg")
    plt.show()


print("Lisco compa")

"""