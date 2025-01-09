import json
import pandas as pd

result_stream = open("", "rb") ### INPUT FILE PARAM (SPECIFY INPUT DIRECTORY)
result_read = result_stream.read()
data_output = json.loads(result_read)
result_stream.close()

hits_sciname = []
hits_accession = []
hits_prot = []
hits_scores = []
hits_identity = []
hits_eval = []
hits_cover = []

hits = data_output["BlastOutput2"][0]["report"]["results"]["search"]["hits"]
for i in range(len(hits)):
    hit = hits[i]
    desc = hit["description"]
    hsps = hit["hsps"][0]
    for j in range(len(desc)):
        organism = desc[j]
        prot = organism["title"].lower()
        sciname = organism["sciname"].lower()
        if "hydrazine dehydrogenase" in prot: ### FILTER NAME PARAM (SPECIFY ENZYME NAME)
            hits_sciname.append(organism["sciname"])
            hits_accession.append(organism["accession"])
            hits_prot.append(organism["title"])

            hits_scores.append(hsps["bit_score"])
            hits_identity.append(hsps["identity"]/hsps["align_len"]*100)
            hits_eval.append(hsps["evalue"])
            cover = abs(hsps["query_from"]-hsps["query_to"]+1)/hsps["align_len"]*100
            hits_cover.append(cover)

sciname_set = list(set(hits_sciname))
sciname_single = []
accession_single = []
prot_single = []
scores_single = []
identity_single = []
eval_single = []
cover_single = []
for sciname in sciname_set:
    idx_sciname = hits_sciname.index(sciname)
    sciname_single.append(hits_sciname[idx_sciname])
    accession_single.append(hits_accession[idx_sciname])
    prot_single.append(hits_prot[idx_sciname])
    scores_single.append(hits_scores[idx_sciname])
    identity_single.append(hits_identity[idx_sciname])
    eval_single.append(hits_eval[idx_sciname])
    cover_single.append(hits_cover[idx_sciname])

dict_filtered = {
    "organism" : sciname_single,
    "accession" : accession_single,
    "prot" : prot_single,
    "score" : scores_single,
    "identity" : identity_single,
    "evalue" : eval_single,
    "cover" : cover_single
}

df_filtered = pd.DataFrame(data=dict_filtered)      
df_filtered.to_csv(".csv") ### OUTPUT FILE PARAM (SPECIFY OUTPUT FILE NAME & DIRECTORY)