import synapseclient
import synapseutils
import pandas as pd
import tarfile
import os

# Login
syn = synapseclient.Synapse()
PROJFOLDER = '~/'
TOKEN = open(PROJFOLDER + 'data/synapse_tokenfile.txt','r').read()
syn.login(authToken=TOKEN)

# Query table
df = pd.read_csv(PROJFOLDER + 'data/protquery_olinkUKB.tsv', sep = '\t')

# Downloading each query
df.apply(
  lambda row: syn.get(
    row['id'],
    downloadLocation = PROJFOLDER + 'data/GWAS_SumStat/' + row['HGNC.symbol']
  ), 
  axis = 1
)

# Untar
df.apply(
  lambda row: tarfile.open(
    PROJFOLDER + 'data/GWAS_SumStat/' + row['HGNC.symbol'] + '/' + row['name']
  ).extractall(
    path = PROJFOLDER + 'data/GWAS_SumStat/' + row['HGNC.symbol']
  ), 
  axis = 1
)

# Remove tar files
df.apply(
  lambda row: os.remove(
    PROJFOLDER + 'data/GWAS_SumStat/' + row['HGNC.symbol'] + '/' + row['name']
  ), 
  axis = 1
)
