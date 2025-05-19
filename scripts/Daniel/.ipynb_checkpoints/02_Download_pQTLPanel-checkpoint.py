import synapseclient
import synapseutils
syn = synapseclient.Synapse()
PROJFOLDER = '/Users/da1078co/Documents/Lund/PhD/Projects/BN_Naeimeh/'
TOKEN = open(PROJFOLDER + 'data/synapse_tokenfile.txt','r').read()
syn.login(authToken=TOKEN)
FOLDER_ID = 'syn51396728'
FILENAME = 'olink_protein_map_3k_v1.tsv'
FILESYN = syn.findEntityId(parent = FOLDER_ID, name = FILENAME)
OUT = PROJFOLDER + 'data/'
syn.get(FILESYN, downloadLocation=OUT)