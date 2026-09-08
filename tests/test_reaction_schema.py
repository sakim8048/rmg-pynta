import pytest
import yaml

from rmgpynta import reaction_schema as rs

VERIFIED_RECORD = {
    "family": "Surface_Dissociation",
    "reactants": [
        {"label": "NX", "adjacency_list": "1 N u0 p1 c0 {2,D} {3,S}\n2 X u0 p0 c0 {1,D}\n3 H u0 p0 c0 {1,S}\n"},
        {"label": "X", "adjacency_list": "1 X u0 p0 c0\n"},
    ],
    "products": [
        {"label": "NHX", "adjacency_list": "1 N u0 p1 c0 {2,T}\n2 X u0 p0 c0 {1,T}\n"},
        {"label": "HX", "adjacency_list": "1 H u0 p0 c0 {2,S}\n2 X u0 p0 c0 {1,S}\n"},
    ],
    "atom_labels_verified": True,
}

UNVERIFIED_RECORD = dict(VERIFIED_RECORD, atom_labels_verified=False)


def test_unmapped_family_raises():
    record = dict(VERIFIED_RECORD, family="Surface_Adsorption_Bidentate")
    with pytest.raises(rs.UnsupportedReactionFamily):
        rs.reaction_record_to_pynta_entry(record, 0)


def test_mapped_family_translates():
    assert rs.translate_family("Surface_Dissociation") == "Dissociation"
    assert rs.translate_family("Surface_Adsorption_Dissociative") == "Dissociative Adsorption"
    assert rs.translate_family("Surface_Abstraction") == "Surface Abstraction"


def test_unverified_labels_raise_by_default():
    with pytest.raises(rs.AtomLabelsNotVerified):
        rs.reaction_record_to_pynta_entry(UNVERIFIED_RECORD, 0)


def test_unverified_labels_allowed_when_overridden():
    entry = rs.reaction_record_to_pynta_entry(UNVERIFIED_RECORD, 0, allow_unverified_labels=True)
    assert entry["reaction_family"] == "Dissociation"


def test_entry_shape_matches_real_pynta_yaml():
    entry = rs.reaction_record_to_pynta_entry(VERIFIED_RECORD, 2)
    assert entry["index"] == 2
    assert entry["reaction"] == "NX + X => NHX + HX"
    assert entry["reaction_family"] == "Dissociation"
    assert entry["reactant"].startswith("multiplicity 1\n")
    assert entry["product"].startswith("multiplicity 1\n")


def test_write_pynta_reactions_yaml_roundtrips(tmp_path):
    out = rs.write_pynta_reactions_yaml([VERIFIED_RECORD, VERIFIED_RECORD], tmp_path / "reaction.yaml")
    loaded = yaml.safe_load(out.read_text())
    assert len(loaded) == 2
    assert loaded[0]["reaction_family"] == "Dissociation"
    assert "reactant" in loaded[0] and "product" in loaded[0]
