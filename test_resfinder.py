import os
import pytest
import pandas as pd
from module.resfinder_result_parser import ResfinderParser, ResfinderCollector

# tests/test_tests.py


@pytest.fixture
def test_data_dir():
    # Adjust this path if your test_data is elsewhere
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "test_data"))

@pytest.mark.parametrize("isolate", ["Isolate01", "Isolate02", "Isolate03", "Isolate04"])
def test_resfinder_parser_basic(test_data_dir, isolate):
    parser = ResfinderParser(test_data_dir, isolate)
    if isolate == "Isolate04":
        assert not parser.has_data
        assert parser.passport.result_summary == "No results found."
        assert parser.resfinder_results_filepath == ""
    else:
        assert parser.has_data
        assert parser.passport.result_summary != ""
        assert os.path.exists(parser.resfinder_results_filepath)
        df = parser.collect_phenotype_results()
        assert isinstance(df, pd.DataFrame)
        assert "isolate_id" in df.columns

@pytest.mark.parametrize("isolate,expected", [
    ("Isolate01", False),
    ("Isolate02", False),
    ("Isolate03", True),
    ("Isolate04", False),
])
def test_resfinder_parser_pointfinder(test_data_dir, isolate, expected):
    parser = ResfinderParser(test_data_dir, isolate)
    assert parser.has_pointfinder_data == expected

def test_resfinder_collector_init(test_data_dir):
    collector = ResfinderCollector(test_data_dir)
    # Should find 4 isolates
    assert set(collector.isolate_dirs) == {"Isolate01", "Isolate02", "Isolate03", "Isolate04"}

def test_resfinder_collector_all_results(test_data_dir):
    collector = ResfinderCollector(test_data_dir)
    results = collector.collect_all_results()
    assert isinstance(results, tuple)
    assert len(results) == 4
    # Each result should be a DataFrame or list/None
    # pointfinder_results, resfinder_results, isolate_summaries, combined_presence_absence
    assert isinstance(results[1], pd.DataFrame)
    assert isinstance(results[2], pd.DataFrame)
    assert isinstance(results[3], pd.DataFrame)
