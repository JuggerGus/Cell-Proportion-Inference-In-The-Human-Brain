# Libraries
import pandas as pd
import numpy as np
import time
import time


    # Read the RNA-seq file. In this work we used RNA-seq from Allen institute. We retrived gene expression from the "donor9861"

# Donor_9861 meta data
head_donor9861=pd.read_csv("SampleAnnot.csv", header=None)
print(head_donor9861)

# Donor_9861 gene expression raw values
conts_donor9861=pd.read_csv("RNAseqCounts.csv", header=None)
print(conts_donor9861)



    # First step, it is only to build the df to the necessary requirements.

col_list = head_donor9861[6].tolist() # Obtain the annotation column of interest
#print(col_list)

conts_donor9861.columns=col_list # Column naming
conts_donor9861=conts_donor9861.transpose()

print("\nRaw matrix:\n")
print(conts_donor9861)

# Save to csv
conts_donor9861.to_csv("donor9861.csv", header=None)





    # Second step, observe the constructed df

sort=pd.read_csv("donor10021.csv", index_col=0)
print(sort)


"""


    # Testing steps



# Testing

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





    # Change raw values to CPM values


# Step 1: Calculate the sum of counts for each sample
sample_totals = np.sum(sort, axis=1)
print(sort)
print(sample_totals)


#print("Librar size:" ,"\n",sample_totals, "\n")

# Step 2: Scale the counts by dividing each count by the total count for the corresponding sample
scaled_matrix = sort / sample_totals[:, np.newaxis]

#print("Scaling factor:" ,"\n",scaled_matrix, "\n")

print(scaled_matrix)
time.sleep(60)

# Step 3: Multiply the scaled counts by 1 million to obtain CPM values
cpm = scaled_matrix * 1e6

print("\nCPM from ChatGPT \n")
cpm=cpm.transpose()
print(cpm)

# Save to csv
cpm.to_csv("donor9861.csv")




    # Third step, sort the df


data=pd.read_csv("CPM_donor9861.csv", index_col=0) 
print(data)

data=data.reindex(sorted(data.columns), axis=1)
print(data)
data.to_csv("CPM_donor9861_sort.csv")



    # Fourth step, convert to .txt

data=pd.read_csv("CPM_donor9861.csv", header=None) 
np.savetxt("CPM_donor9861.txt",data, fmt='%s', delimiter='\t') 





#------------------------------Other way------------------------


# If desired, you can use the bioinfokit.analys library to transform the counts into CPM values.
# In the peper we built the code above to do the normalization step, however these tools
# can be implemented very quickly


"""
# load sugarcane RNA-seq expression dataset (Published in Bedre et al., 2019)
df = get_data('sc_exp').data
print(df.head(2))


# As this data has gene length column, we will drop length column
df = df.drop(['length'], axis=1)
# Make gene column as index column
df = df.set_index('gene')

# Now, normalize raw counts using CPM method 
nm = norm()

# Donor_9861 from Allen Institute
conts_donor9861=pd.read_csv("conts_donor9861_T.csv", index_col=0)
print(conts_donor9861)

nm.cpm(conts_donor9861)

# Get CPM normalized dataframe
cpm_df = nm.cpm_norm
print(cpm_df.head(2))
cpm_df.to_csv("conts_donor9861_Paper.csv")

"""