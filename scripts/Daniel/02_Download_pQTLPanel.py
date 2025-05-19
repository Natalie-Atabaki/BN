import synapseclient
import synapseutils
import pandas as pd

# Login
syn = synapseclient.Synapse()
PROJFOLDER = '/Users/da1078co/Documents/Lund/PhD/Projects/BN_Naeimeh/'
TOKEN = open(PROJFOLDER + 'data/synapse_tokenfile.txt','r').read()
syn.login(authToken=TOKEN)

# Get Olink protein map
FOLDER_ID = 'syn51396728'
FILENAME = 'olink_protein_map_3k_v1.tsv'
FILESYN = syn.findEntityId(parent = FOLDER_ID, name = FILENAME)
OUT = PROJFOLDER + 'data/'
syn.get(FILESYN, downloadLocation=OUT)

# Get file ids
FOLDER_ID = 'syn51365303'
FILEMAP = syn.getChildren(FOLDER_ID, includeTypes=['file'])
FILEDICT = [{'name': item['name'], 'id': item['id']} for item in FILEMAP]
FILEDF = pd.DataFrame(FILEDICT)
FILEDF.to_csv(PROJFOLDER + 'data/olinkfilemapUKB.tsv', sep='\t', index = False)

# Get SNP maps
FOLDER_ID = 'syn51396727'
SNPMAP = syn.getChildren(FOLDER_ID, includeTypes=['file'])
[
  syn.get(
    item['id'], 
    downloadLocation = PROJFOLDER + 'data/olinkukbsnps/'
  ) for item in SNPMAP
]