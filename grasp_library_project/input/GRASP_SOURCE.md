# GRASP sequences (Farley et al., Nucleic Acids Research 2025)

Imported from the deposited GenBank modules (`GRASP_-1.gb` / Assembly Planner).

- Paper: https://academic.oup.com/nar/article/53/20/gkaf1169/8321212
- Data: https://github.com/farleykvdg/GRASP

The CSV junction coordinates record the deposited design. Runtime redesign is
not tied to those indices: every synonymous four-base window fully inside the
invariant ARELF motif (offsets 0–11) is eligible.

Order fragments use configurable physical 5′ and 3′ overhangs at each cloning
level. The deposited defaults are Level −1 `ACAT/ACAA`, Level 0 `CTCA/CTCG`,
and Level 1 `GGAG/AGCG`, with every value written 5′→3′. Retained coding-strand
sites and internal architecture junctions are derived. The deposited workflow
uses pAGM1311 followed by pAGM9121.

Deposited-preset order-strand geometry:
TTTGGTCTCAACAT{{pAGM1311 insert}}TTGTTGAGACCAAA

Deposited PPR block chains (not complete expression constructs):
- 9S: AGGT–CDS1–CTTC–CDS2–TTCG
- 14S: AGGT–CDS1–GTGA–CDS14–CTTC–CDS2–TTCG
- 19S: AGGT–CDS1–GTGA–CDS14–CACG–CDS19–CTTC–CDS2–TTCG
