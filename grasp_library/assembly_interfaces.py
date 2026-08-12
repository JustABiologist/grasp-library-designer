"""Assembly-interface profiles and sequence geometry for GRASP order fragments.

Canonical terminal values name their physical location on the coding-oriented
construct: ``five_prime_end_overhang`` is at its 5' end and
``three_prime_end_overhang`` is at its 3' end. Overhang labels are written
5' to 3', so compatible ends are reverse complements. Where the
bases retained on the assembled coding strand differ from the physical sticky-
end label, a separate ``*_assembled_coding_site`` value makes that explicit.

The custom profile captures requirements only.  It deliberately does not claim
that a vector sequence has been inspected.  The deposited profile captures the
interfaces read from the GRASP planner/deposited constructs, but also stops
short of claiming whole-vector sequence verification.
"""

from __future__ import annotations

import copy
from typing import Any, Mapping, Optional


DNA = frozenset("ACGT")
CANONICAL_NOTATION = "physical_terminal_overhangs_5to3"
FIVE_PRIME_END = "five_prime_end_overhang"
THREE_PRIME_END = "three_prime_end_overhang"
FIVE_PRIME_CODING_SITE = "five_prime_assembled_coding_site"
THREE_PRIME_CODING_SITE = "three_prime_assembled_coding_site"


def reverse_complement(sequence: str) -> str:
    sequence = str(sequence).upper().replace("U", "T")
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _order_fragment_defaults() -> dict[str, Any]:
    return {
        "enzyme": "BsaI",
        "recognition_site": "GGTCTC",
        "clamp_5p": "TTT",
        "spacer_5p": "A",
        "spacer_3p": "T",
        "clamp_3p": "AAA",
    }


def _common_junctions() -> dict[str, dict[str, str]]:
    return {
        "terminal_to_cds2": {
            "upstream_three_prime_end_overhang": "GAAG",
            "downstream_five_prime_end_overhang": "CTTC",
            "assembled_coding_site": "CTTC",
        },
        "cds1_to_cds14": {
            "upstream_three_prime_end_overhang": "TCAC",
            "downstream_five_prime_end_overhang": "GTGA",
            "assembled_coding_site": "GTGA",
        },
        "cds14_to_cds19": {
            "upstream_three_prime_end_overhang": "CGTG",
            "downstream_five_prime_end_overhang": "CACG",
            "assembled_coding_site": "CACG",
        },
    }


def _architectures() -> dict[str, dict[str, list[str]]]:
    return {
        "9S": {
            "blocks": ["CDS1", "CDS2"],
            "joins": ["terminal_to_cds2"],
        },
        "14S": {
            "blocks": ["CDS1", "CDS14", "CDS2"],
            "joins": ["cds1_to_cds14", "terminal_to_cds2"],
        },
        "19S": {
            "blocks": ["CDS1", "CDS14", "CDS19", "CDS2"],
            "joins": [
                "cds1_to_cds14",
                "cds14_to_cds19",
                "terminal_to_cds2",
            ],
        },
    }


def custom_interface_preset() -> dict[str, Any]:
    """Default requirements for a user-supplied vector set.

    No vector sequence is bundled for this profile, so consumers may validate
    the order-fragment geometry and interface requirements, but not the vector
    sequence itself.
    """
    return {
        "profile_name": "custom",
        "notation": CANONICAL_NOTATION,
        "coding_strand_direction": "5prime_to_3prime",
        "order_fragment": _order_fragment_defaults(),
        "level_minus1_entry": {
            "vector_id": "custom_level_minus1_entry",
            FIVE_PRIME_END: "ACAT",
            THREE_PRIME_END: "ACAA",
            FIVE_PRIME_CODING_SITE: "ACAT",
            THREE_PRIME_CODING_SITE: "TTGT",
            "vector_sequence": None,
            "completion_context_5p": None,
            "completion_context_3p": None,
            "release_recognition_site": "GAAGAC",
        },
        "level0": {
            "acceptor_id": "custom_level0_acceptor",
            "acceptor_outer": {
                FIVE_PRIME_END: "CTCA",
                THREE_PRIME_END: "CTCG",
                FIVE_PRIME_CODING_SITE: "CTCA",
                THREE_PRIME_CODING_SITE: "CGAG",
            },
            "ppr_outer": {
                FIVE_PRIME_END: "AGGT",
                THREE_PRIME_END: "CGAA",
                FIVE_PRIME_CODING_SITE: "AGGT",
                THREE_PRIME_CODING_SITE: "TTCG",
            },
            "vector_sequence": None,
        },
        "junctions": _common_junctions(),
        "architectures": _architectures(),
        "final_cassette": {
            "vector_id": "custom_level1_acceptor",
            FIVE_PRIME_END: "GGAG",
            THREE_PRIME_END: "AGCG",
            FIVE_PRIME_CODING_SITE: "GGAG",
            THREE_PRIME_CODING_SITE: "CGCT",
            "vector_sequence": None,
        },
    }


def deposited_grasp_interface_preset() -> dict[str, Any]:
    """Interfaces used by the deposited GRASP assembly planner."""
    return {
        "profile_name": "deposited_grasp",
        "notation": CANONICAL_NOTATION,
        "coding_strand_direction": "5prime_to_3prime",
        "order_fragment": _order_fragment_defaults(),
        "level_minus1_entry": {
            "vector_id": "pAGM1311",
            FIVE_PRIME_END: "ACAT",
            THREE_PRIME_END: "ACAA",
            FIVE_PRIME_CODING_SITE: "ACAT",
            THREE_PRIME_CODING_SITE: "TTGT",
            "vector_sequence": None,
            "completion_context_5p": "GAAG",
            "completion_context_3p": "CTTC",
            "release_recognition_site": "GAAGAC",
        },
        "level0": {
            "acceptor_id": "pAGM9121",
            "acceptor_outer": {
                FIVE_PRIME_END: "CTCA",
                THREE_PRIME_END: "CTCG",
                FIVE_PRIME_CODING_SITE: "CTCA",
                THREE_PRIME_CODING_SITE: "CGAG",
            },
            "ppr_outer": {
                FIVE_PRIME_END: "AGGT",
                THREE_PRIME_END: "CGAA",
                FIVE_PRIME_CODING_SITE: "AGGT",
                THREE_PRIME_CODING_SITE: "TTCG",
            },
            "vector_sequence": None,
        },
        "junctions": _common_junctions(),
        "architectures": _architectures(),
        "final_cassette": {
            "vector_id": "modified_1-1R_pICH47802_lc_p15A_ori_",
            FIVE_PRIME_END: "GGAG",
            THREE_PRIME_END: "AGCG",
            FIVE_PRIME_CODING_SITE: "GGAG",
            THREE_PRIME_CODING_SITE: "CGCT",
            "vector_sequence": None,
        },
    }


PRESETS = {
    "custom": custom_interface_preset,
    "custom_directional_default": custom_interface_preset,
    "deposited": deposited_grasp_interface_preset,
    "deposited_grasp": deposited_grasp_interface_preset,
}


def _deep_update(base: dict[str, Any], update: Mapping[str, Any]) -> None:
    for key, value in update.items():
        if isinstance(value, Mapping) and isinstance(base.get(key), dict):
            _deep_update(base[key], value)
        else:
            base[key] = copy.deepcopy(value)


def _normalize_configured_interfaces(configured: Mapping[str, Any]) -> tuple[str, dict[str, Any]]:
    """Translate the dashboard schema into the canonical internal schema."""
    entry = dict(configured.get("level_minus1_entry", {}))
    level0 = dict(configured.get("level0", {}))
    level1 = dict(configured.get("level1", {}))
    profile_name = str(
        configured.get("preset")
        or entry.get("profile")
        or configured.get("profile_name")
        or "custom"
    )
    if profile_name == "custom_directional_default":
        profile_name = "custom"

    overrides: dict[str, Any] = {
        key: copy.deepcopy(configured[key])
        for key in ("profile_name", "order_fragment", "architectures")
        if key in configured
    }
    source_notation = configured.get("notation", CANONICAL_NOTATION)
    if source_notation != CANONICAL_NOTATION:
        raise ValueError(f"unsupported assembly-interface notation: {source_notation}")
    overrides["notation"] = CANONICAL_NOTATION

    entry_overrides = {
        key: copy.deepcopy(value)
        for key, value in entry.items()
        if key
        in {
            "vector_id",
            FIVE_PRIME_END,
            THREE_PRIME_END,
            FIVE_PRIME_CODING_SITE,
            THREE_PRIME_CODING_SITE,
            "vector_sequence",
            "completion_context_5p",
            "completion_context_3p",
            "release_recognition_site",
        }
    }
    aliases = {
        "vector_name": "vector_id",
    }
    for source, target in aliases.items():
        if source in entry:
            entry_overrides[target] = entry[source]
    if (
        FIVE_PRIME_END in entry_overrides
        and FIVE_PRIME_CODING_SITE not in entry_overrides
    ):
        entry_overrides[FIVE_PRIME_CODING_SITE] = entry_overrides[FIVE_PRIME_END]
    if (
        THREE_PRIME_END in entry_overrides
        and THREE_PRIME_CODING_SITE not in entry_overrides
    ):
        entry_overrides[THREE_PRIME_CODING_SITE] = reverse_complement(
            entry_overrides[THREE_PRIME_END]
        )
    if entry_overrides:
        overrides["level_minus1_entry"] = entry_overrides

    level0_overrides = {
        key: copy.deepcopy(value)
        for key, value in level0.items()
        if key in {"acceptor_id", "acceptor_outer", "ppr_outer", "vector_sequence"}
    }
    for terminal_key in ("acceptor_outer", "ppr_outer"):
        terminal = level0_overrides.get(terminal_key)
        if isinstance(terminal, Mapping):
            terminal = dict(terminal)
            if FIVE_PRIME_END in terminal and FIVE_PRIME_CODING_SITE not in terminal:
                terminal[FIVE_PRIME_CODING_SITE] = terminal[FIVE_PRIME_END]
            if THREE_PRIME_END in terminal and THREE_PRIME_CODING_SITE not in terminal:
                terminal[THREE_PRIME_CODING_SITE] = reverse_complement(
                    terminal[THREE_PRIME_END]
                )
            level0_overrides[terminal_key] = terminal
    if "acceptor_name" in level0:
        level0_overrides["acceptor_id"] = level0["acceptor_name"]
    if level0_overrides:
        overrides["level0"] = level0_overrides

    junction_overrides: dict[str, Any] = copy.deepcopy(
        dict(configured.get("junctions", {}))
    )
    if junction_overrides:
        overrides["junctions"] = junction_overrides

    final_overrides = dict(configured.get("final_cassette", {}))
    if level1:
        for key in (
            FIVE_PRIME_END,
            THREE_PRIME_END,
            FIVE_PRIME_CODING_SITE,
            THREE_PRIME_CODING_SITE,
            "vector_sequence",
        ):
            if key in level1:
                final_overrides[key] = copy.deepcopy(level1[key])
        if "acceptor_name" in level1:
            final_overrides["vector_id"] = level1["acceptor_name"]
    if FIVE_PRIME_END in final_overrides:
        final_overrides[FIVE_PRIME_CODING_SITE] = final_overrides[FIVE_PRIME_END]
    if THREE_PRIME_END in final_overrides:
        final_overrides[THREE_PRIME_CODING_SITE] = reverse_complement(
            final_overrides[THREE_PRIME_END]
        )
    final_overrides = {key: value for key, value in final_overrides.items() if value is not None}
    if final_overrides:
        overrides["final_cassette"] = final_overrides
    return profile_name, overrides


def _validate_dna(value: Any, *, name: str, length: Optional[int] = None) -> str:
    sequence = str(value).upper().replace("U", "T")
    if not sequence or set(sequence) - DNA:
        raise ValueError(f"{name} must contain DNA bases only")
    if length is not None and len(sequence) != length:
        raise ValueError(f"{name} must be {length} nt, got {len(sequence)}")
    return sequence


def validate_assembly_interfaces(profile: Mapping[str, Any]) -> dict[str, Any]:
    """Validate physical sticky ends and their assembled coding-site mapping."""
    result = copy.deepcopy(dict(profile))
    if result.get("notation") != CANONICAL_NOTATION:
        raise ValueError(
            f"assembly interface notation must be {CANONICAL_NOTATION}"
        )

    order = result["order_fragment"]
    for key in ("recognition_site", "clamp_5p", "spacer_5p", "spacer_3p", "clamp_3p"):
        order[key] = _validate_dna(order[key], name=f"order_fragment.{key}")

    entry = result["level_minus1_entry"]
    for key in (
        FIVE_PRIME_END,
        THREE_PRIME_END,
        FIVE_PRIME_CODING_SITE,
        THREE_PRIME_CODING_SITE,
    ):
        entry[key] = _validate_dna(
            entry[key], name=f"level_minus1_entry.{key}", length=4
        )
    if entry[FIVE_PRIME_CODING_SITE] != entry[FIVE_PRIME_END]:
        raise ValueError(
            "level_minus1_entry five-prime assembled coding site must equal "
            "its 5-prime overhang"
        )
    if entry[THREE_PRIME_CODING_SITE] != reverse_complement(entry[THREE_PRIME_END]):
        raise ValueError(
            "level_minus1_entry three-prime assembled coding site must reverse-"
            "complement its 3-prime overhang"
        )
    if entry.get("release_recognition_site") is not None:
        entry["release_recognition_site"] = _validate_dna(
            entry["release_recognition_site"],
            name="level_minus1_entry.release_recognition_site",
        )

    for section_name in ("ppr_outer",):
        section = result["level0"][section_name]
        for key in (
            FIVE_PRIME_END,
            THREE_PRIME_END,
            FIVE_PRIME_CODING_SITE,
            THREE_PRIME_CODING_SITE,
        ):
            section[key] = _validate_dna(
                section[key], name=f"level0.{section_name}.{key}", length=4
            )
        if section[FIVE_PRIME_CODING_SITE] != section[FIVE_PRIME_END]:
            raise ValueError(f"level0.{section_name} 5-prime site mismatch")
        if section[THREE_PRIME_CODING_SITE] != reverse_complement(
            section[THREE_PRIME_END]
        ):
            raise ValueError(f"level0.{section_name} 3-prime site mismatch")
    outer = result["level0"].get("acceptor_outer")
    if outer is not None:
        for key in (
            FIVE_PRIME_END,
            THREE_PRIME_END,
            FIVE_PRIME_CODING_SITE,
            THREE_PRIME_CODING_SITE,
        ):
            outer[key] = _validate_dna(
                outer.get(key), name=f"level0.acceptor_outer.{key}", length=4
            )
        if outer[FIVE_PRIME_CODING_SITE] != outer[FIVE_PRIME_END]:
            raise ValueError("level0.acceptor_outer 5-prime site mismatch")
        if outer[THREE_PRIME_CODING_SITE] != reverse_complement(
            outer[THREE_PRIME_END]
        ):
            raise ValueError("level0.acceptor_outer 3-prime site mismatch")

    final = result["final_cassette"]
    for key in (
        FIVE_PRIME_END,
        THREE_PRIME_END,
        FIVE_PRIME_CODING_SITE,
        THREE_PRIME_CODING_SITE,
    ):
        final[key] = _validate_dna(
            final[key], name=f"final_cassette.{key}", length=4
        )
    if final[FIVE_PRIME_CODING_SITE] != final[FIVE_PRIME_END]:
        raise ValueError("final_cassette 5-prime site mismatch")
    if final[THREE_PRIME_CODING_SITE] != reverse_complement(
        final[THREE_PRIME_END]
    ):
        raise ValueError("final_cassette 3-prime site mismatch")

    junctions = result["junctions"]
    for name, junction in junctions.items():
        upstream_end = _validate_dna(
            junction["upstream_three_prime_end_overhang"],
            name=f"junctions.{name}.upstream_three_prime_end_overhang",
            length=4,
        )
        downstream_end = _validate_dna(
            junction["downstream_five_prime_end_overhang"],
            name=f"junctions.{name}.downstream_five_prime_end_overhang",
            length=4,
        )
        if reverse_complement(upstream_end) != downstream_end:
            raise ValueError(
                f"junctions.{name} physical sticky ends are incompatible: "
                f"reverse_complement({upstream_end}) != {downstream_end}"
            )
        assembled = _validate_dna(
            junction["assembled_coding_site"],
            name=f"junctions.{name}.assembled_coding_site",
            length=4,
        )
        if assembled != downstream_end:
            raise ValueError(
                f"junctions.{name} assembled coding site must equal its downstream "
                "5-prime overhang"
            )
        junction.update(
            upstream_three_prime_end_overhang=upstream_end,
            downstream_five_prime_end_overhang=downstream_end,
            assembled_coding_site=assembled,
        )

    for architecture, layout in result["architectures"].items():
        if len(layout["joins"]) != len(layout["blocks"]) - 1:
            raise ValueError(f"{architecture}: expected one join between each block")
        missing = set(layout["joins"]) - set(junctions)
        if missing:
            raise ValueError(f"{architecture}: undefined junctions {sorted(missing)}")
    return result


def resolve_assembly_interfaces(
    config: Optional[Mapping[str, Any]] = None,
    *,
    preset: Optional[str] = None,
) -> dict[str, Any]:
    """Resolve a named preset plus optional nested configuration overrides."""
    configured: Any = (config or {}).get("assembly_interfaces")
    if isinstance(configured, str):
        preset = configured
        overrides: Mapping[str, Any] = {}
    elif isinstance(configured, Mapping):
        configured_preset, overrides = _normalize_configured_interfaces(configured)
        preset = preset or configured_preset
    elif configured is None:
        overrides = {}
    else:
        raise ValueError("assembly_interfaces must be a preset name or mapping")

    preset = preset or "deposited_grasp"
    try:
        profile = PRESETS[preset]()
    except KeyError as exc:
        raise ValueError(f"unknown assembly-interface preset: {preset}") from exc
    _deep_update(profile, overrides)
    return validate_assembly_interfaces(profile)


def order_fragment_arms(profile: Mapping[str, Any]) -> tuple[str, str]:
    """Return the coding-oriented order strand's BsaI/entry arms.

    The physical 3' sticky-end label is reverse-complemented to obtain the
    bases present at the order strand's 3' assembled site.  The explicit
    assembled-site field is checked against that derivation.
    """
    order = profile["order_fragment"]
    entry = profile["level_minus1_entry"]
    prefix = (
        order["clamp_5p"]
        + order["recognition_site"]
        + order["spacer_5p"]
        + entry[FIVE_PRIME_CODING_SITE]
    )
    suffix = (
        entry[THREE_PRIME_CODING_SITE]
        + order["spacer_3p"]
        + reverse_complement(order["recognition_site"])
        + order["clamp_3p"]
    )
    if entry[FIVE_PRIME_CODING_SITE] != entry[FIVE_PRIME_END]:
        raise ValueError("entry five-prime assembled site does not match its overhang")
    if entry[THREE_PRIME_CODING_SITE] != reverse_complement(entry[THREE_PRIME_END]):
        raise ValueError(
            "entry three-prime assembled site must reverse-complement its overhang"
        )
    return prefix, suffix


def build_order_fragment(payload: str, profile: Mapping[str, Any]) -> str:
    payload = _validate_dna(payload, name="payload")
    prefix, suffix = order_fragment_arms(profile)
    return prefix + payload + suffix


def extract_order_payload(sequence: str, profile: Mapping[str, Any]) -> str:
    sequence = _validate_dna(sequence, name="order fragment")
    prefix, suffix = order_fragment_arms(profile)
    if not sequence.startswith(prefix) or not sequence.endswith(suffix):
        raise ValueError("order fragment does not match the configured BsaI/entry arms")
    payload = sequence[len(prefix) : -len(suffix)]
    if len(payload) < 8:
        raise ValueError("released payload is too short to contain two interfaces")
    return payload
