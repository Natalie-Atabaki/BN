import synapseclient
import synapseutils
syn = synapseclient.Synapse()
TOKEN="ABC123" #This should be changed accordingly
syn.login(authToken=TOKEN)
FOLDER_ID = "syn51824537"

TRAIT="GLCGN"
FILENAME = "Glucagon_4891_50.txt.gz"
FILESYN = syn.findEntityId(parent=FOLDER_ID, name = FILENAME)
OUT="/Users/da1078co/Documents/PhD/Projects/BN_Naeimeh/data/GWAS_SumStat/" + TRAIT
syn.get(FILESYN, downloadLocation=OUT)

TRAIT="GLP1R"
FILENAME = "GLP1R_13085_18.txt.gz"
FILESYN = syn.findEntityId(parent=FOLDER_ID, name = FILENAME)
OUT="/Users/da1078co/Documents/PhD/Projects/BN_Naeimeh/data/GWAS_SumStat/" + TRAIT
syn.get(FILESYN, downloadLocation=OUT)