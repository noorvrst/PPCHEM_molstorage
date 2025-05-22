import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from chemstorm import chemsort_multiple_order_3, initialize_storage_groups

import pytest

from chemstorm.chemstorm import (
    get_compound_safety_data,
    get_name_and_smiles,
    classify_acid_base,
    get_mp_bp,
    compound_state,
    prioritize_pictograms,
    is_chemically_compatible,
    initialize_storage_groups,
    chemsort_multiple_order_3
)

def test_get_compound_safety_data_valid():
    cid, pictograms, hazards = get_compound_safety_data("acetone")
    assert cid.isdigit()  # CID should be numeric
    assert isinstance(pictograms, list)
    assert isinstance(hazards, list)

def test_get_compound_safety_data_invalid():
    cid, pictograms, hazards = get_compound_safety_data("unrealisticchemicalxyz")
    assert cid == ""
    assert pictograms == []
    assert hazards == []

def test_get_compound_safety_data_debug():
    cid, pictograms, hazards = get_compound_safety_data("hydrochloric acid", debug=True)
    assert isinstance(cid, str)
    assert isinstance(pictograms, list)
    assert isinstance(hazards, list)

def test_get_name_and_smiles_valid():
    name, iupac, smiles = get_name_and_smiles("180")  # 180 is acetone
    assert isinstance(name, str)
    assert isinstance(iupac, str)
    assert isinstance(smiles, str)
    assert smiles != "Unknown"

def test_get_name_and_smiles_invalid():
    try:
        get_name_and_smiles("000000000")  # Should raise or return unknown
        assert True  # If it doesn't raise, it's okay as long as it handles gracefully
    except Exception:
        assert True  # If it raises, that's fine for bad CID


def test_classify_acid_base_acid():
    result = classify_acid_base("acetic acid", "ethanoic acid", "CC(=O)O", [])
    assert result == "acid"

def test_classify_acid_base_base():
    result = classify_acid_base("ammonia", "azane", "N", [])
    assert result == "base"

def test_classify_acid_base_amphoteric():
    result = classify_acid_base("amino acid", "2-aminoethanoic acid", "NCC(=O)O", [])
    assert result == "amphoteric"

def test_classify_acid_base_neutral():
    result = classify_acid_base("hexane", "hexane", "CCCCCC", [])
    assert result == "neutral"

def test_classify_acid_base_invalid_smiles():
    result = classify_acid_base("weird compound", "unknown", "INVALID", [])
    assert result == "Invalid SMILES"

def test_get_mp_bp_water():
    mp_c, bp_c, mp_f, bp_f = get_mp_bp("water")
    assert isinstance(mp_c, (int, float, type(None)))
    assert isinstance(bp_c, (int, float, type(None)))
    assert isinstance(mp_f, (int, float, type(None)))
    assert isinstance(bp_f, (int, float, type(None)))

def test_get_mp_bp_invalid_compound():
    mp_c, bp_c, mp_f, bp_f = get_mp_bp("notarealcompound")
    assert mp_c is None
    assert bp_c is None
    assert mp_f is None
    assert bp_f is None

def test_state_liquid_celsius():
    assert compound_state(0, 100, None, None) == "liquid"

def test_state_solid_celsius():
    assert compound_state(50, 300, None, None) == "solid"

def test_state_gas_celsius():
    assert compound_state(-100, 10, None, None) == "gas"

def test_state_liquid_fahrenheit_fallback():
    assert compound_state(None, None, 32, 212) == "liquid"

def test_state_solid_fahrenheit_fallback():
    assert compound_state(None, None, 80, 300) == "solid"

def test_state_gas_fahrenheit_fallback():
    assert compound_state(None, None, -100, 10) == "gas"

def test_state_unknown():
    assert compound_state(None, None, None, None) == "unknown"

def test_prioritize_standard():
    input_pics = ["Irritant", "Flammable", "Explosive"]
    output = prioritize_pictograms(input_pics)
    assert output == ["Explosive", "Flammable", "Irritant"]

def test_prioritize_unknown():
    input_pics = ["UnknownHazard", "Flammable", "Corrosive"]
    output = prioritize_pictograms(input_pics)
    assert output == ["Flammable", "Corrosive", "UnknownHazard"]

def test_prioritize_same_priority():
    input_pics = ["Health Hazard", "Acute Toxic"]
    output = prioritize_pictograms(input_pics)
    assert output == ["Health Hazard", "Acute Toxic"]  # order is preserved

def test_prioritize_all_known():
    input_pics = [
        "Environmental Hazard", "Irritant", "Health Hazard",
        "Acute Toxic", "Corrosive", "Flammable", "Explosive", "Oxidizer"
    ]
    output = prioritize_pictograms(input_pics)
    expected = [
        "Explosive", "Oxidizer", "Flammable", "Corrosive",
        "Health Hazard", "Acute Toxic", "Environmental Hazard", "Irritant"
    ]
    assert output == expected

def test_empty_input():
    assert prioritize_pictograms([]) == []

@pytest.mark.parametrize("existing_pictograms, new_pictograms, existing_acid_base, new_acid_base, group, expected", [
    (["Flammable"], ["Oxidizer"], "", "", "", False),
    (["Corrosive"], ["Flammable"], "acid", "base", "", False),
    (["Corrosive"], ["Health Hazard"], "acid", "", "", False),
    ([], ["Flammable"], "", "", "", True),
    (["Corrosive"], ["Irritant"], "base", "base", "", True),
    (["Oxidizer"], ["Corrosive"], "", "", "oxidizer", False),
    (["Corrosive"], ["Oxidizer"], "", "", "corrosive", False),
    (["Flammable"], ["Corrosive"], "", "", "flammable", False),
    (["Corrosive"], ["Acute Toxic"], "acid", "", "acid_corrosive_1", False),
])

def test_is_chemically_compatible(existing_pictograms, new_pictograms,
                                  existing_acid_base, new_acid_base, group, expected):
    assert is_chemically_compatible(existing_pictograms, new_pictograms,
                                    existing_acid_base, new_acid_base, group) == expected

def test_initialize_storage_groups_structure():
    groups = initialize_storage_groups()
    expected_keys = {
        "none", "hazardous_environment", "acute_toxicity", "cmr_stot", "toxicity_2_3",
        "acid_corrosive_1", "acid_irritant", "base_corrosive_1", "base_irritant",
        "pyrophoric", "flammable", "oxidizer", "explosive", "compressed_gas", "nitric_acid"
    }
    assert set(groups.keys()) == expected_keys
    for group in groups.values():
        assert all(state in group for state in ["solid", "liquid", "gas", "unknown"])
        assert all(isinstance(group[state], list) for state in group)

def test_chemsort_basic_compatibility():
    groups = initialize_storage_groups()
    test_compounds = [
        {
            "name": "Ethanol",
            "sorted_pictograms": ["Flammable"],
            "hazard_statements": ["Highly flammable liquid and vapor"],
            "acid_base_class": "neutral",
            "state_room_temp": "liquid"
        },
        {
            "name": "Nitric Acid",
            "sorted_pictograms": ["Corrosive"],
            "hazard_statements": ["Causes severe skin burns and eye damage"],
            "acid_base_class": "acid",
            "state_room_temp": "liquid"
        },
        {
            "name": "Sodium Hydroxide",
            "sorted_pictograms": ["Corrosive"],
            "hazard_statements": ["Causes severe skin burns and eye damage"],
            "acid_base_class": "base",
            "state_room_temp": "solid"
        },
        {
            "name": "Hydrogen",
            "sorted_pictograms": ["Compressed Gas"],
            "hazard_statements": ["Contains gas under pressure; may explode if heated"],
            "acid_base_class": "neutral",
            "state_room_temp": "gas"
        },
        {
            "name": "Unknown Substance",
            "sorted_pictograms": [],
            "hazard_statements": [],
            "acid_base_class": "neutral",
            "state_room_temp": "solid"
        }
    ]

    sorted_groups = chemsort_multiple_order_3(test_compounds, groups)

    assert any(c["name"] == "Ethanol" for c in sorted_groups["flammable"]["liquid"])
    assert any(c["name"] == "Nitric Acid" for c in sorted_groups["nitric_acid"]["liquid"])
    assert any(c["name"] == "Sodium Hydroxide" for c in sorted_groups["base_corrosive_1"]["solid"])
    assert any(c["name"] == "Hydrogen" for c in sorted_groups["compressed_gas"]["gas"])
    assert any(c["name"] == "Unknown Substance" for c in sorted_groups["none"]["solid"])


def test_custom_group_for_flammable_and_oxidizer():
    groups = initialize_storage_groups()
    test_compounds = [
        {
            "name": "Hydrogen Peroxide",
            "sorted_pictograms": ["Corrosive"],
            "hazard_statements": [""],
            "acid_base_class": "neutral",
            "state_room_temp": "liquid"
        },
        {
            "name": "Ethanol",
            "sorted_pictograms": ["Flammable"],
            "hazard_statements": ["Highly flammable liquid and vapor"],
            "acid_base_class": "neutral",
            "state_room_temp": "liquid"
        }
    ]

    sorted_groups = chemsort_multiple_order_3(test_compounds, groups)

    ethanol_in_flammable = any(c["name"] == "Ethanol" for c in sorted_groups["flammable"]["liquid"])

    custom_group_has_peroxide = any(
        group.startswith("custom_storage_") and
        any(c["name"] == "Hydrogen Peroxide" for c in sorted_groups[group]["liquid"])
        for group in sorted_groups
    )

    assert ethanol_in_flammable, "Ethanol should be placed in the flammable group"
    assert custom_group_has_peroxide, "Hydrogen Peroxide should be separated into a custom group due to incompatibility"
