import pandas as pd
import pytest

from grasp_library.binder import annotate_module_roles, describe_part_id
from grasp_library.import_grasp import compile_target_gap


@pytest.mark.parametrize(
    ("part_id", "role", "upstream", "downstream", "fifth", "last", "rna"),
    [
        ("B_LD5N", "B", "A", "C", "N", "D", "U"),
        ("C_LN5T", "C", "B", "D", "T", "N", "A"),
        ("D_LD5T", "D", "C", "E", "T", "D", "G"),
        ("1A_5T_AGGT", "A", "", "B", "T", "", "A or G"),
        ("1A_5N_AATG", "A", "", "B", "N", "", "C or U"),
        ("2A_LN5N", "A", "", "B", "N", "N", "C"),
        ("1E_LD5N", "E", "D", "", "N", "D", "U"),
        ("2E_LD", "E", "D", "", "", "D", "G or U"),
        ("14A_LN5T", "A", "", "B", "T", "N", "A"),
        ("19E_LN5N", "E", "D", "", "N", "N", "C"),
    ],
)
def test_describe_part_id_slot_and_ppr_code(
    part_id, role, upstream, downstream, fifth, last, rna
):
    info = describe_part_id(part_id)
    assert info["assembly_role"] == role
    assert info["joins_upstream_role"] == upstream
    assert info["joins_downstream_role"] == downstream
    assert info["ppr_5th_aa"] == fifth
    assert info["ppr_last_aa"] == last
    assert info["target_rna_base"] == rna


def test_annotate_puts_role_columns_beside_part_id():
    table = pd.DataFrame(
        {
            "optimized_part_id": ["B_LD5N_v1"],
            "part_id": ["B_LD5N"],
            "oh5": ["ACTC"],
        }
    )
    out = annotate_module_roles(table)
    assert list(out.columns[:8]) == [
        "optimized_part_id",
        "part_id",
        "assembly_role",
        "joins_upstream_role",
        "joins_downstream_role",
        "ppr_5th_aa",
        "ppr_last_aa",
        "target_rna_base",
    ]
    assert out.loc[0, "assembly_role"] == "B"
    assert out.loc[0, "target_rna_base"] == "U"


def test_assembly_plan_labels_five_part_roles():
    plan = compile_target_gap("UUACACGUG", architecture="9S")
    assert plan["assembly_role"].tolist() == ["A", "B", "C", "D", "E"] * 2
    cds1 = plan[plan["assembly_group"] == "CDS1"]
    assert cds1["joins_downstream_role"].tolist() == ["B", "C", "D", "E", ""]
    assert "U" in set(plan["target_rna_base"])
