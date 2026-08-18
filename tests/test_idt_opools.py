from grasp_library.idt_opools import price_idt_opool


def test_opool_floor_covers_typical_oneshot_pool():
    quote = price_idt_opool([114] * 8 + [84, 90, 91])
    assert quote["n_oligos"] == 11
    assert quote["total_bases"] < 3300
    assert quote["dna_eur"] == 109.00
    assert quote["total_eur"] == 109.00
    assert quote["eligible_scales"] == ["10 pmol", "50 pmol"]
    assert quote["phospho_eur"] == 17.93


def test_opool_second_tier_starts_after_3300_bases():
    quote = price_idt_opool([3300, 1])
    assert quote["total_bases"] == 3301
    assert quote["dna_eur"] == 109.04


def test_opool_phospho_adds_per_oligo():
    quote = price_idt_opool([80, 80], phosphorylate_5prime=True)
    assert quote["dna_eur"] == 109.00
    assert quote["phospho_eur"] == 3.26
    assert quote["total_eur"] == 112.26
