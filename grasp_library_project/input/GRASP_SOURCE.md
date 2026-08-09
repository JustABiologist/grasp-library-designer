# GRASP sequences (Farley et al., Nucleic Acids Research 2025)

Imported from the deposited GenBank modules (`GRASP_-1.gb` / Assembly Planner).

- Paper: https://academic.oup.com/nar/article/53/20/gkaf1169/8321212
- Data: https://github.com/farleykvdg/GRASP

Cut indices (`junction_map.mask_start_0based`) are taken from the GenBank
overhang features relative to each part’s in-frame coding window. Protein
sequences are translations of those windows (synonymous redesign only).

Default 9S overhang set: AGGT–ACTC–AAGA–GCAC–TGAA–CTTC–ACTC–AAGA–GCAC–TGAA–TTCG
(1A AGGT N-terminal fusion; AATG variants included as alternate parts).
