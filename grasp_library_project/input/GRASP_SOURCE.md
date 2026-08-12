# GRASP sequences (Farley et al., Nucleic Acids Research 2025)

Imported from the deposited GenBank modules (`GRASP_-1.gb` / Assembly Planner).

- Paper: https://academic.oup.com/nar/article/53/20/gkaf1169/8321212
- Data: https://github.com/farleykvdg/GRASP

The GenBank overhang features define the native design. Redesign candidates are
not tied to those indices: the program enumerates every synonymous four-base
window whose motif-relative start is 0–11 inside invariant `ARELF`, then
rematerializes the part window and coding mask (synonymous redesign only).

Order fragments use configurable BsaI entry-vector overhangs; BpiI then releases
modules into the configured Level 0 acceptor. The deposited GRASP preset uses
pAGM1311 then pAGM9121. Terminal-side labels name the physical end of the
coding-oriented construct: N-terminal side = 5′ end and C-terminal side = 3′
end. Every overhang label is written 5′→3′, so compatible ends are reverse
complements. Custom defaults use entry 5′/N—3′/C `AACA/GGAG`, CDS1
3′/C—CDS2 5′/N `CTTC/GAAG`, and final cassette 5′/N—3′/C `GCCC/GCGA`.

Deposited pAGM1311 preset order-strand geometry:
TTTGGTCTCAACAT{{pAGM1311 insert}}TTGTTGAGACCAAA

Deposited GRASP PPR block chains (not complete expression cassettes):
- 9S: AGGT–CDS1–CTTC–CDS2–TTCG
- 14S: AGGT–CDS1–GTGA–CDS14–CTTC–CDS2–TTCG
- 19S: AGGT–CDS1–GTGA–CDS14–CACG–CDS19–CTTC–CDS2–TTCG
