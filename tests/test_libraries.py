from rmgpynta.libraries import update_input_py_libraries

TEMPLATE = """
database(
    thermoLibraries=[
        'surfaceThermoPt111',
        'primaryThermoLibrary',
        'DFT_QCI_thermo',
    ],

    reactionLibraries=[],
    seedMechanisms=[],
    kineticsFamilies=['surface', 'default'],
)
"""


def test_inserts_new_names_and_keeps_existing(tmp_path):
    template_path = tmp_path / "input.py"
    template_path.write_text(TEMPLATE)
    out_path = tmp_path / "round_001" / "input.py"

    update_input_py_libraries(
        template_path, out_path, thermo_library_names=["round_000"], kinetics_library_names=["round_000"]
    )
    text = out_path.read_text()

    assert "'round_000'" in text
    assert "'surfaceThermoPt111'" in text
    assert "'primaryThermoLibrary'" in text
    assert "'DFT_QCI_thermo'" in text
    # kineticsFamilies untouched
    assert "kineticsFamilies=['surface', 'default']" in text


def test_does_not_duplicate_already_listed_name(tmp_path):
    template_path = tmp_path / "input.py"
    template_path.write_text(TEMPLATE)
    out_path = tmp_path / "input.py"

    update_input_py_libraries(
        template_path, out_path, thermo_library_names=["surfaceThermoPt111"], kinetics_library_names=[]
    )
    text = out_path.read_text()
    assert text.count("'surfaceThermoPt111'") == 1
