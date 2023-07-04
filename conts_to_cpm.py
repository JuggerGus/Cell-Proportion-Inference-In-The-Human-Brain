# Librerias 
import pandas as pd
import numpy as np
import time


# Paths
path_donor9861="/Users/gusano2398gmail.com/Desktop/Tesis/Data/Data_Allen/Bulk/rnaseq_donor9861/"
path_donor10021="/Users/gusano2398gmail.com/Desktop/Tesis/Data/Data_Allen/Bulk/rnaseq_donor10021/"

# Donor_9861
head_donor9861=pd.read_csv(path_donor9861+"SampleAnnot.csv", header=None)
conts_donor9861=pd.read_csv(path_donor9861+"RNAseqCounts.csv", header=None)

# Donor_10021
conts_donor10021=pd.read_csv(path_donor10021+"RNAseqCounts.csv", header=None)
head_donor10021=pd.read_csv(path_donor10021+"SampleAnnot.csv", header=None)





# Primer paso
"""
col_list = head_donor10021[6].tolist() #Obtener la columna de anotaciones de interes #CAMBIAR
#print(col_list)
conts_donor10021.columns=col_list #Colocar nombre de columnas #CAMBIAR
conts_donor10021=conts_donor10021.transpose() #CAMBIAR
print("\nRaw matrix:\n")
print(conts_donor10021) #CAMBIAR
# Guardar en Csv
conts_donor10021.to_csv("donor10021.csv", header=None) #CAMBIAR


"""



# Segundo paso


sort=pd.read_csv("donor10021.csv", index_col=0) #CAMBIAR
print(sort)



"""

# Comprobar

# Example gene expression matrix
expression_matrix = np.array([[10, 10, 100],
                              [10, 20, 200],
                              [10, 30, 300]])

matrix=pd.DataFrame(expression_matrix)



m=[]

for sample in sort.index:
    values=sort.loc[sample, sort.columns]
    #print("Values:", sample,"\n",values, "\n")

    library_sizes = np.sum(values)
    #print("Librar size:", sample ,"\n",library_sizes, "\n")
    scaling_factor = library_sizes / 1e6
   #print("Scaling factor:",sample ,"\n",scaling_factor, "\n")

    for v in values:
        cpm_matrix = v / scaling_factor
        #print(cpm_matrix)
        #time.sleep(5)
        m.append(cpm_matrix)
        #time.sleep(5)

ma=pd.DataFrame(m)
print("CPM matrix\n")
print(ma)

"""





"""
# Step 1: Calculate the sum of counts for each sample
sample_totals = np.sum(sort, axis=1)
#print("Librar size:" ,"\n",sample_totals, "\n")
# Step 2: Scale the counts by dividing each count by the total count for the corresponding sample
scaled_matrix = sort / sample_totals[:, np.newaxis]
#print("Scaling factor:" ,"\n",scaled_matrix, "\n")
# Step 3: Multiply the scaled counts by 1 million to obtain CPM values
cpm = scaled_matrix * 1e6

print("\nCPM from ChatGPT \n")
cpm=cpm.transpose()
print(cpm)

cpm.to_csv("CPM_donor10021.csv") #CAMBIAR



"""



# Tercer paso
"""
data=pd.read_csv("CPM_donor10021.csv", index_col=0) #CAMBIAR
print(data)
data=data.reindex(sorted(data.columns), axis=1)
print(data)
data.to_csv("CPM_donor10021_sort.csv")

"""

# Cuarto Paso

data=pd.read_csv("CPM_donor10021_sort.csv", header=None) #CAMBIAR
np.savetxt("CPM_donor10021.txt",data, fmt='%s', delimiter='\t') #CAMBIAR
