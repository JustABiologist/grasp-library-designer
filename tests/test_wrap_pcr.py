from Bio.SeqUtils import MeltingTemp as melting

from grasp_library.dna import reverse_complement
from grasp_library.gga_split import wrap_geometry, wrap_payload, wrap_pcr_primers


def test_bsai_wrap_arms_are_pool_pcr_primers_at_55c():
    geometry = wrap_geometry("BsaI")
    primers = wrap_pcr_primers(geometry)
    oligo = wrap_payload("AATG" + "ACGT" * 10 + "AAGC", geometry)

    assert primers["forward"] == "ACAGCCAGGTCTCA"
    assert primers["reverse"] == "TATCGGCGGTCTCA"
    assert primers["forward"] != primers["reverse"]
    assert primers["tm_forward_c"] >= 55.0
    assert primers["tm_reverse_c"] >= 55.0
    assert oligo.startswith(primers["forward"])
    assert oligo.endswith(reverse_complement(primers["reverse"]))
    assert oligo.count("GGTCTC") == 1
    observed = melting.Tm_NN(
        primers["forward"],
        nn_table=melting.DNA_NN3,
        Na=50,
        Mg=1.5,
        dnac1=250,
        dnac2=0,
        dNTPs=0.2,
    )
    assert observed >= 55.0
