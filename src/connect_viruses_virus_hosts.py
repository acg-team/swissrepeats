import csv
import re

f_sp = "/home/matteo/polybox/MSc_ACLS/swissrepeat/data/swissprot_annotations.tsv"
f_tr = "/home/matteo/polybox/MSc_ACLS/swissrepeat/results/tr_annotations/concatenated.csv"

# 1. Extract sp ids with trs.

with open(f_tr, "r") as fh:
    csv_reader = csv.reader(fh, delimiter=',', quotechar='|')
    tr_header = next(csv_reader)
    sp_id = set([row[0] for row in csv_reader])
    # Extract species:


# 2. Extract species & virus hosts from f_sp
with open(f_sp, "r") as fh:
    csv_reader = csv.reader(fh, delimiter='\t', quotechar='|')
    sp_header = next(csv_reader)
    sp_annotations = [row for row in csv_reader]

virus_host_index = sp_header.index("Virus hosts")
entry_index = sp_header.index("Entry")
species_id_index = sp_header.index("Taxonomic lineage IDs (SPECIES)")
species_name_index = sp_header.index("Taxonomic lineage (SPECIES)")

# 2.a) Extract virus hosts Tax IDs and Tax names from f_sp
virus_hosts = [[host for host in row[virus_host_index].split(";") if host != ""] for row in sp_annotations]

p_tax_id = re.compile("TaxID: (\d+)")
virus_hosts_TaxID = [[p_tax_id.search(host).group(1) for host in row] if len(row) != 0 else [] for row in virus_hosts]
virus_hosts_TaxID_flat = [item for sublist in virus_hosts_TaxID for item in sublist]

p_tax_name = re.compile("([\w\s()]+)")
virus_hosts_TaxName = [[p_tax_name.search(host).group(1).strip() for host in row] if len(row) != 0 else [] for row in virus_hosts]
virus_hosts_TaxName_flat = [item for sublist in virus_hosts_TaxName for item in sublist]

# Map these to Swiss-Prot proteins
virus_host_TaxID_proteins = {i: [] for i in set(virus_hosts_TaxID_flat)}
for row, hosts in zip(sp_annotations, virus_hosts_TaxID):
    for host in hosts:
        virus_host_TaxID_proteins[host].append(row[entry_index])

virus_host_TaxName_proteins = {i: [] for i in set(virus_hosts_TaxName_flat)}
for row, hosts in zip(sp_annotations, virus_hosts_TaxName):
    for host in hosts:
        virus_host_TaxName_proteins[host].append(row[entry_index])



# 2.b) species Tax IDs and Tax names from f_sp
species_TaxID = [row[species_id_index] for row in sp_annotations]
species_TaxName = [row[species_name_index] for row in sp_annotations]
species_TaxID_to_TaxName = {id:name for id, name in zip(species_TaxID, species_TaxName)}
species_TaxName_to_TaxID = {name:id for id, name in zip(species_TaxID, species_TaxName)}


# Map these to Swiss-Prot proteins
species_TaxID_proteins = {i: [] for i in set(species_TaxID)}
for row in sp_annotations:
    species_TaxID_proteins[row[species_id_index]].append(row[entry_index])

species_TaxName_proteins = {i: [] for i in set(species_TaxName)}
for row in sp_annotations:
    species_TaxName_proteins[row[species_name_index]].append(row[entry_index])


# Check for overlaps
virus_hosts_in_swissprot_by_TaxID = set(virus_hosts_TaxID_flat) & set(species_TaxID)
# -> 73 overlaps for the TaxID

virus_hosts_in_swissprot_by_TaxName = set(virus_hosts_TaxName_flat) & set(species_TaxName)
# -> 339 overlaps for the TaxName

# Total overlap
virus_hosts_in_swissprot_total = set([species_TaxID_to_TaxName[i] for i in virus_hosts_in_swissprot_by_TaxID]) | virus_hosts_in_swissprot_by_TaxName
# -> 342
virus_hosts_in_swissprot_diff_TaxName = set([species_TaxID_to_TaxName[i] for i in virus_hosts_in_swissprot_by_TaxID]) - virus_hosts_in_swissprot_by_TaxName
virus_hosts_in_swissprot_diff_TaxID = [species_TaxName_to_TaxID[i] for i in virus_hosts_in_swissprot_diff_TaxName]

# {'Odocoileus virginianus (White-tailed deer)', 'Varecia variegata (Black-and-white ruffed lemur) (Lemur variegatus)', "Saccharomyces cerevisiae (Baker's yeast)"}


# 3. Check which species have most annotations.

virus_hosts_in_swissprot_by_TaxName_sorted = list(virus_hosts_in_swissprot_by_TaxName)
virus_hosts_in_swissprot_by_TaxName_sorted.sort(key=lambda x: -len(virus_host_TaxName_proteins[x]))

for vh in virus_hosts_in_swissprot_by_TaxName_sorted:
    print("{}: {}".format(vh, len(virus_host_TaxName_proteins[vh])))

hosts_well_represented_by_virus_sequences = {vh: len(proteins) for vh, proteins in virus_host_TaxName_proteins.items() if len(proteins) >= 100}
hosts_well_represented_by_own_sequence = {vh: len(virus_host_TaxName_proteins[vh]) for vh in virus_hosts_in_swissprot_by_TaxName_sorted if len(virus_host_TaxName_proteins[vh]) >= 100}

# How good is he overlap between hosts_well_represented_by_virus_sequences and hosts_well_represented_by_own_sequence?
len(set(hosts_well_represented_by_virus_sequences.keys()) & set(hosts_well_represented_by_own_sequence.keys()))


for vh in virus_hosts_in_swissprot_diff_TaxID:
    print("{}: {}".format(species_TaxID_to_TaxName[vh], len(virus_host_TaxID_proteins[vh])))


print("\n".join(sorted(hosts_well_represented_by_own_sequence.keys())))


# 4. Print all pairs of viruses and hosts
virus_host_pairs = set()
for sp_virushosts, sp_virus_TaxName, sp_virus_TaxID in zip(virus_hosts_TaxName, species_TaxName, species_TaxID):
    for sp_virushost_TaxName in sp_virushosts:
        virus_host_pairs.add((sp_virus_TaxName, sp_virus_TaxID, sp_virushost_TaxName))


virus_host_pairs_host_taxname_in_sp = [i for i in virus_host_pairs if i[2] in virus_hosts_in_swissprot_by_TaxName]

results_file = "/home/matteo/polybox/MSc_ACLS/swissrepeat/results/virus_host_pairs_host_taxname_in_sp.csv"
with open(results_file, "w") as fh:
    fh.write(",".join(["virus_TaxName", "virus_TaxID", "virushost_TaxName"]))
    for i in virus_host_pairs_host_taxname_in_sp:
        fh.write("\n"+ ",".join(i))