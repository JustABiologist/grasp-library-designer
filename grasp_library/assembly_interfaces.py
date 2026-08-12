"""Assembly-interface profiles and sequence geometry for GRASP order fragments.

All terminal interfaces use ``directional_terminal_5p`` notation: the N value
is the 5' sequence on the assembled coding strand and the C value is the 5'
sequence on the opposite strand.  Therefore adjacent terminals are compatible
when ``reverse_complement(upstream_C) == downstream_N``.

The custom profile captures requirements only.  It deliberately does not claim
that a vector sequence has been inspected.  The deposited profile captures the
interfaces read from the GRASP planner/deposited constructs, but also stops
short of claiming whole-vector sequence verification.
"""

from __future__ import annotations

import copy
from typing import Any, Mapping, Optional


DNA = frozenset("ACGT")


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
            "upstream_c_5p": "CTTC",
            "downstream_n_5p": "GAAG",
            "assembled_plus_site": "CTTC",
        },
        "cds1_to_cds14": {
            "upstream_c_5p": "GTGA",
            "downstream_n_5p": "TCAC",
            "assembled_plus_site": "GTGA",
        },
        "cds14_to_cds19": {
            "upstream_c_5p": "CACG",
            "downstream_n_5p": "CGTG",
            "assembled_plus_site": "CACG",
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
        "profile_name": "custom_directional_default",
        "notation": "directional_terminal_5p",
        "order_fragment": _order_fragment_defaults(),
        "level_minus1_entry": {
            "vector_id": "custom_level_minus1_entry",
            "n_overhang_5p": "AACA",
            "c_overhang_5p": "GGAG",
            "vector_sequence": None,
            "completion_context_5p": None,
            "completion_context_3p": None,
            "release_recognition_site": "GAAGAC",
        },
        "level0": {
            "acceptor_id": "custom_level0_acceptor",
            "acceptor_outer": None,
            "ppr_outer": {
                "n_overhang_5p": "AGGT",
                "c_overhang_5p": "TTCG",
            },
            "vector_sequence": None,
        },
        "junctions": _common_junctions(),
        "architectures": _architectures(),
        "final_cassette": {
            "vector_id": "custom_level1_acceptor",
            "n_overhang_5p": "GCCC",
            "c_overhang_5p": "GCGA",
            "vector_sequence": None,
        },
    }


def deposited_grasp_interface_preset() -> dict[str, Any]:
    """Interfaces used by the deposited GRASP assembly planner."""
    return {
        "profile_name": "deposited_grasp",
        "notation": "directional_terminal_5p",
        "order_fragment": _order_fragment_defaults(),
        "level_minus1_entry": {
            "vector_id": "pAGM1311",
            "n_overhang_5p": "ACAT",
            "c_overhang_5p": "ACAA",
            "vector_sequence": None,
            "completion_context_5p": "GAAG",
            "completion_context_3p": "CTTC",
            "release_recognition_site": "GAAGAC",
        },
        "level0": {
            "acceptor_id": "pAGM9121",
            "acceptor_outer": {
                "n_overhang_5p": "CTCA",
                "c_overhang_5p": "CGAG",
            },
            "ppr_outer": {
                "n_overhang_5p": "AGGT",
                "c_overhang_5p": "TTCG",
            },
            "vector_sequence": None,
        },
        "junctions": _common_junctions(),
        "architectures": _architectures(),
        "final_cassette": {
            "vector_id": "modified_1-1R_pICH47802_lc_p15A_ori_",
            "n_overhang_5p": "GGAG",
            "c_overhang_5p": "CGCT",
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
        or "custom_directional_default"
    )
    if profile_name == "custom":
        profile_name = "custom_directional_default"

    overrides: dict[str, Any] = {
        key: copy.deepcopy(configured[key])
        for key in ("profile_name", "order_fragment", "architectures")
        if key in configured
    }
    notation = configured.get("notation", configured.get("overhang_notation"))
    if notation is not None:
        overrides["notation"] = notation

    entry_overrides = {
        key: value
        for key, value in entry.items()
        if key in {"vector_id", "n_overhang_5p", "c_overhang_5p", "vector_sequence", "completion_context_5p", "completion_context_3p", "release_recognition_site"}
    }
    aliases = {
        "vector_name": "vector_id",
        "n_terminal_overhang": "n_overhang_5p",
        "c_terminal_overhang": "c_overhang_5p",
    }
    for source, target in aliases.items():
        if source in entry:
            entry_overrides[target] = entry[source]
    if entry_overrides:
        overrides["level_minus1_entry"] = entry_overrides

    level0_overrides = {
        key: value
        for key, value in level0.items()
        if key in {"acceptor_id", "acceptor_outer", "ppr_outer", "vector_sequence"}
    }
    if "acceptor_name" in level0:
        level0_overrides["acceptor_id"] = level0["acceptor_name"]
    if "acceptor_n_terminal_overhang" in level0 or "acceptor_c_terminal_overhang" in level0:
        default_outer = level0_overrides.get("acceptor_outer") or {}
        level0_overrides["acceptor_outer"] = {
            "n_overhang_5p": level0.get(
                "acceptor_n_terminal_overhang", default_outer.get("n_overhang_5p")
            ),
            "c_overhang_5p": level0.get(
                "acceptor_c_terminal_overhang", default_outer.get("c_overhang_5p")
            ),
        }
    if level0_overrides:
        overrides["level0"] = level0_overrides

    junction_overrides: dict[str, Any] = copy.deepcopy(
        dict(configured.get("junctions", {}))
    )
    for legacy_name, values in level0.get("block_junctions", {}).items():
        canonical_name = (
            "terminal_to_cds2" if legacy_name == "cds1_to_cds2" else legacy_name
        )
        existing = dict(values)
        junction_overrides[canonical_name] = {
            "upstream_c_5p": existing.get(
                "upstream_c_5p", existing.get("upstream_c")
            ),
            "downstream_n_5p": existing.get(
                "downstream_n_5p", existing.get("downstream_n")
            ),
            "assembled_plus_site": existing.get(
                "assembled_plus_site",
                existing.get("upstream_c_5p", existing.get("upstream_c")),
            ),
        }
    if junction_overrides:
        overrides["junctions"] = junction_overrides

    final_overrides = dict(configured.get("final_cassette", {}))
    if level1:
        final_overrides.update(
            {
                "vector_id": level1.get("acceptor_name", final_overrides.get("vector_id")),
                "n_overhang_5p": level1.get(
                    "n_terminal_overhang", final_overrides.get("n_overhang_5p")
                ),
                "c_overhang_5p": level1.get(
                    "c_terminal_overhang", final_overrides.get("c_overhang_5p")
                ),
            }
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
    """Return a normalized profile after enforcing directional-pair geometry."""
    result = copy.deepcopy(dict(profile))
    if result.get("notation") != "directional_terminal_5p":
        raise ValueError("assembly interface notation must be directional_terminal_5p")

    order = result["order_fragment"]
    for key in ("recognition_site", "clamp_5p", "spacer_5p", "spacer_3p", "clamp_3p"):
        order[key] = _validate_dna(order[key], name=f"order_fragment.{key}")

    entry = result["level_minus1_entry"]
    for key in ("n_overhang_5p", "c_overhang_5p"):
        entry[key] = _validate_dna(entry[key], name=f"level_minus1_entry.{key}", length=4)
    if entry.get("release_recognition_site") is not None:
        entry["release_recognition_site"] = _validate_dna(
            entry["release_recognition_site"],
            name="level_minus1_entry.release_recognition_site",
        )

    for section_name in ("ppr_outer",):
        section = result["level0"][section_name]
        for key in ("n_overhang_5p", "c_overhang_5p"):
            section[key] = _validate_dna(section[key], name=f"level0.{section_name}.{key}", length=4)
    outer = result["level0"].get("acceptor_outer")
    if outer is not None:
        for key in ("n_overhang_5p", "c_overhang_5p"):
            outer[key] = _validate_dna(
                outer.get(key), name=f"level0.acceptor_outer.{key}", length=4
            )

    final = result["final_cassette"]
    for key in ("n_overhang_5p", "c_overhang_5p"):
        final[key] = _validate_dna(final[key], name=f"final_cassette.{key}", length=4)

    junctions = result["junctions"]
    for name, junction in junctions.items():
        upstream = _validate_dna(junction["upstream_c_5p"], name=f"junctions.{name}.upstream_c_5p", length=4)
        downstream = _validate_dna(junction["downstream_n_5p"], name=f"junctions.{name}.downstream_n_5p", length=4)
        if reverse_complement(upstream) != downstream:
            raise ValueError(
                f"junctions.{name} is not a directional terminal pair: "
                f"reverse_complement({upstream}) != {downstream}"
            )
        assembled = _validate_dna(junction.get("assembled_plus_site", upstream), name=f"junctions.{name}.assembled_plus_site", length=4)
        junction.update(
            upstream_c_5p=upstream,
            downstream_n_5p=downstream,
            assembled_plus_site=assembled,
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

    preset = preset or "custom_directional_default"
    try:
        profile = PRESETS[preset]()
    except KeyError as exc:
        raise ValueError(f"unknown assembly-interface preset: {preset}") from exc
    _deep_update(profile, overrides)
    return validate_assembly_interfaces(profile)


def order_fragment_arms(profile: Mapping[str, Any]) -> tuple[str, str]:
    order = profile["order_fragment"]
    entry = profile["level_minus1_entry"]
    prefix = (
        order["clamp_5p"]
        + order["recognition_site"]
        + order["spacer_5p"]
        + entry["n_overhang_5p"]
    )
    suffix = (
        reverse_complement(entry["c_overhang_5p"])
        + order["spacer_3p"]
        + reverse_complement(order["recognition_site"])
        + order["clamp_3p"]
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
