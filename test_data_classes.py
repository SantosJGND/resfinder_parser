import pytest
import pandas as pd
from module.data_classes import Phenotype, IsolatePhenotypes, SeqRegion, IsolateSummary


@pytest.fixture
def empty_isolate_phenotypes():
    return IsolatePhenotypes("test_isolate", "No resistance found")


@pytest.fixture
def single_phenotype_isolate():
    iso = IsolatePhenotypes("test_isolate", "Resistance detected")
    phenotype = Phenotype("ampicillin")
    phenotype.amr_classes = "beta_lactam"
    phenotype.amr_resistant = "True"
    phenotype.amr_species_relevant = "True"
    phenotype.grade = "1"
    phenotype.seq_regions = []
    iso.add_phenotype(phenotype)
    return iso


@pytest.fixture
def multi_phenotype_isolate():
    iso = IsolatePhenotypes("test_isolate", "Multiple resistances")

    p1 = Phenotype("ampicillin")
    p1.amr_classes = "beta_lactam"
    p1.amr_resistant = "True"
    p1.amr_species_relevant = "True"
    p1.grade = "1"
    p1.seq_regions = []

    p2 = Phenotype("ciprofloxacin")
    p2.amr_classes = "fluoroquinolone"
    p2.amr_resistant = "True"
    p2.amr_species_relevant = "True"
    p2.grade = "2"
    p2.seq_regions = []

    iso.add_phenotype(p1)
    iso.add_phenotype(p2)
    return iso


@pytest.fixture
def phenotype_with_list_fields():
    iso = IsolatePhenotypes("test_isolate", "Resistance detected")
    phenotype = Phenotype("ampicillin")
    phenotype.amr_classes = ["beta_lactam", "penicillin"]
    phenotype.amr_resistant = "True"
    phenotype.amr_species_relevant = "True"
    phenotype.grade = "1"
    phenotype.seq_regions = []
    iso.add_phenotype(phenotype)
    return iso


@pytest.fixture
def phenotype_with_empty_list():
    iso = IsolatePhenotypes("test_isolate", "Resistance detected")
    phenotype = Phenotype("ampicillin")
    phenotype.amr_classes = []
    phenotype.amr_resistant = "True"
    phenotype.amr_species_relevant = "True"
    phenotype.grade = "1"
    phenotype.seq_regions = []
    iso.add_phenotype(phenotype)
    return iso


@pytest.fixture
def phenotype_with_none_field():
    iso = IsolatePhenotypes("test_isolate", "Resistance detected")
    phenotype = Phenotype("ampicillin")
    phenotype.amr_classes = None
    phenotype.amr_resistant = "True"
    phenotype.amr_species_relevant = "True"
    phenotype.grade = "1"
    phenotype.seq_regions = []
    iso.add_phenotype(phenotype)
    return iso


@pytest.fixture
def phenotype_with_seq_regions():
    iso = IsolatePhenotypes("test_isolate", "Resistance detected")
    phenotype = Phenotype("ampicillin")
    phenotype.amr_classes = "beta_lactam"
    phenotype.amr_resistant = "True"
    phenotype.amr_species_relevant = "True"
    phenotype.grade = "1"
    phenotype.seq_regions = [
        SeqRegion("blaTEM-1;;1;;AF456215"),
        SeqRegion("blaOXA-1;;1;;X04560")
    ]
    iso.add_phenotype(phenotype)
    return iso


@pytest.fixture
def phenotype_with_multiple_seq_regions_multiple_phenotypes():
    iso = IsolatePhenotypes("test_isolate", "Multiple resistances")

    p1 = Phenotype("ampicillin")
    p1.amr_classes = "beta_lactam"
    p1.amr_resistant = "True"
    p1.amr_species_relevant = "True"
    p1.grade = "1"
    p1.seq_regions = [
        SeqRegion("blaTEM-1;;1;;AF456215"),
        SeqRegion("blaOXA-1;;1;;X04560")
    ]

    p2 = Phenotype("ciprofloxacin")
    p2.amr_classes = "fluoroquinolone"
    p2.amr_resistant = "True"
    p2.amr_species_relevant = "True"
    p2.grade = "2"
    p2.seq_regions = [
        SeqRegion("gyrA;;1;;NC_000001")
    ]

    iso.add_phenotype(p1)
    iso.add_phenotype(p2)
    return iso


class TestPhenotypeDataFrame:
    def test_empty_phenotypes_returns_dataframe_with_columns(self, empty_isolate_phenotypes):
        df = empty_isolate_phenotypes.phenotype_dataframe()
        
        assert isinstance(df, pd.DataFrame)
        expected_columns = ["antibiotic", "amr_classes", "amr_resistant", 
                           "amr_species_relevant", "grade", "seq_regions"]
        assert list(df.columns) == expected_columns

    def test_empty_phenotypes_returns_empty_dataframe(self, empty_isolate_phenotypes):
        df = empty_isolate_phenotypes.phenotype_dataframe()
        
        assert df.empty

    def test_single_phenotype_returns_one_row(self, single_phenotype_isolate):
        df = single_phenotype_isolate.phenotype_dataframe()
        
        assert len(df) == 1
        assert df.iloc[0]["antibiotic"] == "ampicillin"
        assert df.iloc[0]["amr_classes"] == "beta_lactam"
        assert df.iloc[0]["amr_resistant"] == "True"

    def test_multiple_phenotypes_returns_multiple_rows(self, multi_phenotype_isolate):
        df = multi_phenotype_isolate.phenotype_dataframe()
        
        assert len(df) == 2
        antibiotics = df["antibiotic"].tolist()
        assert "ampicillin" in antibiotics
        assert "ciprofloxacin" in antibiotics

    def test_list_field_non_empty_joined_with_semicolon(self, phenotype_with_list_fields):
        df = phenotype_with_list_fields.phenotype_dataframe()
        
        assert df.iloc[0]["amr_classes"] == "beta_lactam;penicillin"

    def test_list_field_empty_returns_zero_string(self, phenotype_with_empty_list):
        df = phenotype_with_empty_list.phenotype_dataframe()
        
        assert df.iloc[0]["amr_classes"] == "0"

    def test_none_field_returns_none_string(self, phenotype_with_none_field):
        df = phenotype_with_none_field.phenotype_dataframe()
        
        assert df.iloc[0]["amr_classes"] == "None"

    def test_seq_regions_preserves_seqregion_objects(self, phenotype_with_seq_regions):
        df = phenotype_with_seq_regions.phenotype_dataframe()
        
        seq_regions_value = df.iloc[0]["seq_regions"]
        assert isinstance(seq_regions_value, str)
        assert "blaTEM-1" in seq_regions_value
        assert "blaOXA-1" in seq_regions_value

    def test_phenotype_dataframe_columns_order(self, single_phenotype_isolate):
        df = single_phenotype_isolate.phenotype_dataframe()
        
        expected_columns = ["antibiotic", "amr_classes", "amr_resistant", 
                           "amr_species_relevant", "grade", "seq_regions"]
        assert list(df.columns) == expected_columns


class TestAllGenes:
    def test_empty_isolate_returns_empty_list(self, empty_isolate_phenotypes):
        genes = empty_isolate_phenotypes.all_genes()
        
        assert genes == []

    def test_single_gene_returns_list_with_gene(self, phenotype_with_seq_regions):
        genes = phenotype_with_seq_regions.all_genes()
        
        assert len(genes) == 2
        assert "blaTEM-1" in genes
        assert "blaOXA-1" in genes

    def test_multiple_phenotypes_multiple_genes(self, phenotype_with_multiple_seq_regions_multiple_phenotypes):
        genes = phenotype_with_multiple_seq_regions_multiple_phenotypes.all_genes()
        
        assert len(genes) == 3
        assert "blaTEM-1" in genes
        assert "blaOXA-1" in genes
        assert "gyrA" in genes


class TestCollectAllGenesAffected:
    def test_empty_isolate_returns_empty_dict(self, empty_isolate_phenotypes):
        empty_isolate_phenotypes.collect_all_genes_affected()
        
        assert empty_isolate_phenotypes.all_genes_affected == {}

    def test_single_phenotype_single_gene(self, phenotype_with_seq_regions):
        phenotype_with_seq_regions.collect_all_genes_affected()
        
        result = phenotype_with_seq_regions.all_genes_affected
        assert "blaTEM-1" in result
        assert "blaOXA-1" in result
        assert result["blaTEM-1"] == ["ampicillin"]
        assert result["blaOXA-1"] == ["ampicillin"]

    def test_gene_associated_with_multiple_phenotypes(self, phenotype_with_multiple_seq_regions_multiple_phenotypes):
        phenotype_with_multiple_seq_regions_multiple_phenotypes.collect_all_genes_affected()
        
        result = phenotype_with_multiple_seq_regions_multiple_phenotypes.all_genes_affected
        assert "blaTEM-1" in result
        assert "blaOXA-1" in result
        assert "gyrA" in result


class TestGeneAffected:
    def test_gene_not_present_returns_empty_string(self, single_phenotype_isolate):
        result = single_phenotype_isolate.gene_affected("nonexistent_gene")
        
        assert result == ""

    def test_gene_present_returns_phenotype_names(self, phenotype_with_seq_regions):
        phenotype_with_seq_regions.collect_all_genes_affected()
        
        result = phenotype_with_seq_regions.gene_affected("blaTEM-1")
        
        assert result == "ampicillin"

    def test_gene_associated_with_multiple_phenotypes(self, phenotype_with_multiple_seq_regions_multiple_phenotypes):
        phenotype_with_multiple_seq_regions_multiple_phenotypes.collect_all_genes_affected()
        
        result = phenotype_with_multiple_seq_regions_multiple_phenotypes.gene_affected("blaTEM-1")
        
        assert "ampicillin" in result


class TestSeqRegion:
    def test_parsing_gene_version_homolog(self):
        sr = SeqRegion("blaTEM-1;;1;;AF456215")
        
        assert sr.gene == "blaTEM-1"
        assert sr.version == "1"
        assert sr.homolog == "AF456215"

    def test_string_representation(self):
        sr = SeqRegion("blaTEM-1;;1;;AF456215")
        
        assert str(sr) == "blaTEM-1;;1;;AF456215"

    def test_repr_equals_str(self):
        sr = SeqRegion("blaTEM-1;;1;;AF456215")
        
        assert repr(sr) == str(sr)


class TestIsolateSummary:
    def test_dataclass_creation(self):
        summary = IsolateSummary(
            key="version_key",
            provided_species="Escherichia coli",
            result_summary="Resistance detected",
            databases="ResFinder-4.1.0"
        )
        
        assert summary.key == "version_key"
        assert summary.provided_species == "Escherichia coli"
        assert summary.result_summary == "Resistance detected"
        assert summary.databases == "ResFinder-4.1.0"
