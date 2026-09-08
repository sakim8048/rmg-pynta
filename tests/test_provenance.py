from pathlib import Path

from rmgpynta import provenance

FIXTURE = Path(__file__).parent / "fixtures" / "sample_chem_annotated_surface.inp"


def test_parses_all_five_reactions():
    records = provenance.parse_annotated_chemkin(FIXTURE)
    assert len(records) == 5


def test_rate_rule_average_node_classified():
    records = provenance.parse_annotated_chemkin(FIXTURE)
    r = records[0]  # H(8)+NH2(D)(9)<=>NH3(2), "Estimated from node ... in family ..."
    assert r["kinetics_source"] == "rate_rule"
    assert r["rate_rule_path"] == "Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing"
    assert r["family"] == "R_Recombination"
    assert r["chemkin_index"] == 1
    assert r["rmg_index"] == 3
    assert r["flux_pairs"] == ["NH2(D)(9), NH3(2)", "H(8), NH3(2)"]


def test_training_reaction_classified():
    records = provenance.parse_annotated_chemkin(FIXTURE)
    r = records[1]  # H(8)+H(8)<=>H2(3), "Matched reaction 56 ..."
    assert r["kinetics_source"] == "training_reaction"
    assert r["training_reaction"]["number"] == 56
    assert r["training_reaction"]["family"] == "R_Recombination"
    assert r["rate_rule_path"] == "Root_1R->H_N-2R->S_2CHNO->H"


def test_surface_reaction_with_stick_trailer():
    records = provenance.parse_annotated_chemkin(FIXTURE)
    r = records[2]  # X(5)+CH4(1)<=>CH4X(13), STICK trailer
    assert r["kinetics_source"] == "training_reaction"
    assert r["family"] == "Surface_Adsorption_vdW"
    assert r["trailer"] == ["STICK"]
    assert r["reaction_string"] == "X(5)+CH4(1)<=>CH4X(13)"


def test_rate_rule_average_with_degeneracy():
    records = provenance.parse_annotated_chemkin(FIXTURE)
    r = records[3]  # dissociative adsorption, "Average of [...]" + degeneracy 3.0
    assert r["kinetics_source"] == "rate_rule_average"
    assert r["rate_rule_path"] == "From training reaction 1 used for H2;VacantSite1;VacantSite2"
    assert r["degeneracy"] == 3.0


def test_library_reaction_classified():
    records = provenance.parse_annotated_chemkin(FIXTURE)
    r = records[4]  # "! Library reaction: pt111_bootstrap_round_000"
    assert r["kinetics_source"] == "library"
    assert r["library_name"] == "pt111_bootstrap_round_000"


def test_summarize_rollup():
    records = provenance.parse_annotated_chemkin(FIXTURE)
    summary = provenance.summarize(records)
    assert summary["n_reactions"] == 5
    assert summary["by_source"]["library"] == 1
    assert summary["by_source"]["training_reaction"] == 2
    assert summary["fraction_library_backed"] == 1 / 5
