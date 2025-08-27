from dataclasses import dataclass
from typing import List, Dict
import pandas as pd


@dataclass
class IsolateSummary:
    key: str
    provided_species: str
    result_summary: str
    databases: str


class SeqRegion:
    def __init__(self, seq_region_output: str):
        seq_region_output = seq_region_output.split(";;")

        self.gene = seq_region_output[0].strip()
        self.version = seq_region_output[1].strip()
        self.homolog = seq_region_output[2].strip()

    def __str__(self):
        return f"{self.gene};;{self.version};;{self.homolog}"

    def __repr__(self):
        return str(self)


class Phenotype:
    phenotype_fields = [
        "amr_classes",
        "amr_resistant",
        "amr_species_relevant",
        "grade",
        "seq_regions",
    ]

    def __init__(self, name):
        self.name = name
        self.seq_regions: List[SeqRegion] = []

    def add_seq_region(self, seq_region: str):
        self.seq_regions.append(SeqRegion(seq_region))

    def add_regions(self, regions: List[str]):
        for region in regions:
            self.add_seq_region(region)

    def feed_data(self, data: dict):
        for field in self.phenotype_fields:
            found = data.get(field, None)
            if found is not None:
                if isinstance(found, list):
                    if field == "seq_regions":
                        self.add_regions(found)
                        continue
                    found = ";".join(found)

            setattr(self, field, str(found))


class IsolatePhenotypes:
    def __init__(self, isolate_id: str, result_summary: str):
        self.isolate_id = isolate_id
        self.result_summary = result_summary
        self.phenotypes: Dict[str, Phenotype] = {}
        self.all_genes_affected: Dict[str, list] = {}

    def add_phenotype(self, phenotype: Phenotype):
        self.phenotypes[phenotype.name] = phenotype

    def phenotype_dataframe(self):
        data_keep = []
        for phenotype_name, phenotype in self.phenotypes.items():
            phenotype_data = [phenotype_name]
            for field in phenotype.phenotype_fields:
                found = getattr(phenotype, field, None)
                if found is not None:
                    if isinstance(found, list):
                        if len(found) == 0:
                            found = "0"
                        else:
                            found = ";".join([str(x) for x in found])

                phenotype_data.append(str(found))

            data_keep.append(phenotype_data)

        df = pd.DataFrame(
            data_keep, columns=["antibiotic"] + phenotype.phenotype_fields
        )

        return df

    def collect_all_genes_affected(self):
        for phenotype_name, phenotype in self.phenotypes.items():
            for seq_region in phenotype.seq_regions:
                if seq_region.gene not in self.all_genes_affected.keys():
                    self.all_genes_affected[seq_region.gene] = [phenotype_name]
                else:
                    if phenotype_name not in self.all_genes_affected[seq_region.gene]:
                        self.all_genes_affected[seq_region.gene].append(phenotype_name)

    def all_genes(self):
        genes = []
        for _, phenotype in self.phenotypes.items():
            for seq_region in phenotype.seq_regions:
                genes.append(seq_region.gene)

        return genes

    def gene_affected(self, gene: str):
        gene_present = self.all_genes_affected.get(gene, None)

        if gene_present is None:
            return ""
        else:
            return "; ".join(gene_present)
