"""D16: the two DN2 selector fixes (reference-proteome query; species-level match).  No network."""
import select_dn2 as SEL

HITS = [["UP000077096", "Fervidobacterium ngatamarikiense", "1", "2000"],
        ["UP000002524", "Deinococcus radiodurans (strain ATCC 13939 / DSM 20539 / R1)", "243230", "3085"],
        ["UP000001425", "Synechocystis sp. (strain ATCC 27184 / PCC 6803 / Kazusa)", "1111708", "3507"],
        ["UP000099999", "Synechocystis sp. PCC 6714", "2", "3400"]]


def test_query_uses_reference_true():
    q = SEL.proteome_query("Deinococcus radiodurans")
    assert "reference:true" in q and "proteome_type" not in q
    assert q == '(organism_name:"Deinococcus radiodurans") AND (reference:true)'


def test_species_level_match():
    m = SEL.match_reference_proteome("Deinococcus radiodurans", HITS)
    assert m["upid"] == "UP000002524" and m["n_proteins"] == 3085


def test_genus_only_has_no_reference_proteome():
    # bug: "Fervidobacterium" alone used to match another species' proteome
    assert SEL.match_reference_proteome("Fervidobacterium", HITS) is None
    assert SEL.genus_species("Fervidobacterium") is None


def test_placeholder_epithet_is_genus_only():
    # DEVIATIONS item 13 (stricter): "sp." is not a species
    assert SEL.match_reference_proteome("Synechocystis sp. pcc 6803", HITS) is None
    assert SEL.genus_species("Synechocystis sp.") is None
    assert SEL.genus_species("Bacillus spp") is None


def test_no_hit_is_none():
    assert SEL.match_reference_proteome("Frankia alni", HITS) is None
