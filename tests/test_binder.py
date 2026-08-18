from grasp_library.binder import describe_binder, rna_to_binder_aa


def test_library_tract_starts_at_the_solvating_helix():
    aa = rna_to_binder_aa("UUACACGUG")
    assert aa.startswith("QGGNSEEPRKSFDERPERGVVS")
    assert not aa.startswith("M")


def test_oneshot_orf_prepends_the_start_codon():
    aa = rna_to_binder_aa("UUACACGUG", include_start_codon=True)
    info = describe_binder("UUACACGUG")
    assert aa.startswith("MQGGNSEEPRKSFDERPERGVVS")
    assert info["aa_sequence"] == aa
    assert info["include_start_codon"] is True
    assert info["aa_length"] == len(rna_to_binder_aa("UUACACGUG")) + 1
