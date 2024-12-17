from typing import Optional
import pandas as pd
import json
import os
import datetime
import logging
from dataclasses import dataclass
from typing import List, Dict
from module.data_classes import IsolatePhenotypes, Phenotype, IsolateSummary


class ResfinderParser:
    pointfinder_results_filename: str = "PointFinder_results.txt"
    resfinder_results_filename = "ResFinder_results_tab.txt"
    resfinder_results_suffix = ".json"

    def __init__(self, RESFINDER_dir, isolate_dir):
        self.isolate_id = isolate_dir

        self.pointfinder_results_filepath = os.path.join(
            RESFINDER_dir, isolate_dir, self.pointfinder_results_filename
        )
        self.resfinder_results_filepath = os.path.join(
            RESFINDER_dir,
            isolate_dir,
            f"{self.isolate_id}{self.resfinder_results_suffix}",
        )

        self.time_now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        self.has_data = os.path.isfile(
            self.pointfinder_results_filepath
        ) and os.path.isfile(self.resfinder_results_filepath)

        if self.has_data:
            self.passport = self.collect_summary()
        else:
            self.passport = self.empty_passport

        self.phenotypes = IsolatePhenotypes(
            self.isolate_id, self.passport.result_summary
        )

    def resfinder_json_parse_antibiotics(self):
        """Parse the resfinder json output into a pandas dataframe."""
        content = self.read_json()
        ## info to get:

        for phenotype_str, data in content["phenotypes"].items():
            antibiotic = Phenotype(phenotype_str)
            antibiotic.feed_data(data)
            self.phenotypes.add_phenotype(antibiotic)

    @staticmethod
    def resfinder_json_summary(content) -> IsolateSummary:
        """
        Parse the resfinder json output into a IsolateSummary object."""
        analysis_key = content["key"]
        provided_species = content["provided_species"]
        result_summary = content["result_summary"]
        databases: dict = content["databases"]
        databases_to_keep = ["resfinder", "pointfinder"]
        databases = "; ".join(
            [
                f"{v['database_name']}-{v['database_version']}"
                for k, v in databases.items()
                if v["database_name"].lower() in databases_to_keep
            ]
        )
        return IsolateSummary(analysis_key, provided_species, result_summary, databases)

    @property
    def empty_passport(self):
        return IsolateSummary("", "", "No results found.", "")

    def read_json(self):
        with open(self.resfinder_results_filepath, "r") as f:
            results = json.load(f)

        return results

    def collect_summary(self) -> IsolateSummary:
        results = self.read_json()
        return self.resfinder_json_summary(results)

    def add_analyis_columns(self, df: pd.DataFrame):
        df["isolate_id"] = self.isolate_id
        df["analysis_date"] = self.time_now
        df["version"] = self.passport.key

        # place the isolate_id, analysis_date and resfinder version columns at the front
        cols = df.columns.tolist()
        cols = cols[-3:] + cols[:-3]
        df = df[cols]

        return df

    def collect_pointfinder_results(self):
        results = pd.read_csv(self.pointfinder_results_filepath, sep="\t")

        results = self.add_analyis_columns(results)

        return results

    def collect_resfinder_results(self):
        self.resfinder_json_parse_antibiotics()
        isolate_results = self.phenotypes.phenotype_dataframe()
        isolate_results = self.add_analyis_columns(isolate_results)

        return isolate_results

    def seq_regions_parse(self):
        results = self.read_json()

        seq_reg = []
        for reg, map in results["seq_regions"].items():
            for phenotype in map["phenotypes"]:
                seq_reg.append([reg, phenotype, map["query_id"], map["identity"]])

        seq_reg_df = pd.DataFrame(
            seq_reg, columns=["seq_region", "antibiotic", "query_id", "identity"]
        )

        return seq_reg_df

    def extend_resfinder_results(self, resfinder_df):

        seq_reg_df = self.seq_regions_parse()

        def contig_id_info(phenotype: pd.DataFrame):
            pheno_df = seq_reg_df[seq_reg_df.antibiotic == phenotype]

            if pheno_df.shape[0] == 0:
                return pd.Series(["", 0])

            output = pd.Series(
                ["; ".join(pheno_df.query_id.values), pheno_df.identity.values[0]]
            )

            return output

        resfinder_df[["contigs", "identity"]] = resfinder_df.antibiotic.apply(
            contig_id_info
        )

        return resfinder_df

    def isolate_summary(self):
        summary_df = pd.DataFrame(
            [
                [
                    self.passport.databases,
                    self.passport.provided_species,
                    self.passport.result_summary,
                ]
            ],
            columns=["databases", "provided_species", "result_summary"],
        )

        summary_df = self.add_analyis_columns(summary_df)

        return summary_df

    @property
    def empty_isolate_summary(self):
        return pd.DataFrame(
            [
                [
                    self.isolate_id,
                    self.time_now,
                    self.passport.key,
                    self.passport.databases,
                    "",
                    "No results found.",
                ]
            ],
            columns=[
                "isolate_id",
                "analysis_date",
                "version",
                "databases",
                "provided_species",
                "result_summary",
            ],
        )


class ResfinderCollector:
    def __init__(self, RESFINDER_dir: str, output_dir: Optional[str] = None):
        self.RESFINDER_dir = RESFINDER_dir
        self.output_dir = output_dir
        self.isolate_dirs = self.resfinder_dirs()
        self.logger = logging.getLogger(__name__)
        # set level to info to get the info messages
        self.logger.setLevel(logging.INFO)
        # log to console
        self.logger.addHandler(logging.StreamHandler())
        self.logger.info(
            f"Found {len(self.isolate_dirs)} isolate directories in {self.RESFINDER_dir}"
        )

    def _is_isolate_dir(self, dirpath: str):
        if dirpath.startswith("."):
            return False

        if os.path.isdir(os.path.join(self.RESFINDER_dir, dirpath)) is False:
            return False

        return True

    def resfinder_dirs(self):
        isolate_dirs = os.listdir(self.RESFINDER_dir)
        isolate_dirs = [x for x in isolate_dirs if self._is_isolate_dir(x)]
        return isolate_dirs

    def genes_affected(self, isolate_phenotypes: List[IsolatePhenotypes]):
        all_genes = [
            gene for isolate in isolate_phenotypes for gene in isolate.all_genes()
        ]
        all_genes = list(set(all_genes))

        data_to_keep = []

        for isolate in isolate_phenotypes:
            isolate_line = [isolate.isolate_id, isolate.result_summary]
            isolate.collect_all_genes_affected()
            for gene in all_genes:
                isolate_line.append(isolate.gene_affected(gene))

            data_to_keep.append(isolate_line)

        df = pd.DataFrame(
            data_to_keep, columns=["isolate_id", "result_summary"] + all_genes
        )

        return df

    def pointfinder_summary(
        self, pointfinder_results: pd.DataFrame
    ) -> Optional[pd.DataFrame]:

        if pointfinder_results.empty:
            return None

        isolate_ids = pointfinder_results["isolate_id"].unique()

        pointfinder_known = pointfinder_results[
            pointfinder_results["Resistance"] != "Unknown"
        ]

        groups = []
        # for _, group in pointfinder_results.groupby("isolate_id"):
        for isolate_id in isolate_ids:
            group = pointfinder_known[pointfinder_known["isolate_id"] == isolate_id]

            if group.empty is False:
                ## merge rows with the same mutation, but different resistance. concatenate the resistance
                group = group.groupby(["isolate_id", "Mutation"]).agg(
                    {"Resistance": lambda x: "; ".join(x)}
                )
                group = group.reset_index()

                matrix_df = group.pivot(
                    index="isolate_id", columns="Mutation", values="Resistance"
                )

                matrix_df = matrix_df.fillna("")
                matrix_df = matrix_df.reset_index()

                groups.append(matrix_df)
            else:
                groups.append(
                    pd.DataFrame(
                        [[isolate_id] + [""] * (len(group.columns) - 1)],
                        columns=group.columns,
                    )
                )

        matrix_df = pd.concat(groups, axis=0)

        return matrix_df

    def collect_all_results(self):
        pointfinder_results = []
        resfinder_results = []
        isolate_summaries = []
        isolate_phenotypes = []

        for isolate_dir in os.listdir(self.RESFINDER_dir):
            if isolate_dir.startswith("."):
                continue

            if os.path.isdir(os.path.join(self.RESFINDER_dir, isolate_dir)) is False:
                continue

            isolate_parser = ResfinderParser(self.RESFINDER_dir, isolate_dir)

            if isolate_parser.has_data is False:
                self.logger.info(
                    f"No results found for isolate {isolate_parser.isolate_id}"
                )
                isolate_summary = isolate_parser.empty_isolate_summary
                isolate_summaries.append(isolate_summary)

                continue

            isolate_pointfinder_results = isolate_parser.collect_pointfinder_results()
            isolate_resfinder_results = isolate_parser.collect_resfinder_results()
            isolate_resfinder_results = isolate_parser.extend_resfinder_results(
                isolate_resfinder_results
            )
            isolate_summary = isolate_parser.isolate_summary()

            pointfinder_results.append(isolate_pointfinder_results)
            resfinder_results.append(isolate_resfinder_results)
            isolate_summaries.append(isolate_summary)
            isolate_phenotypes.append(isolate_parser.phenotypes)

        if len(isolate_phenotypes) == 0:
            return None, None, None, None

        pointfinder_results = pd.concat(pointfinder_results, axis=0)
        resfinder_results = pd.concat(resfinder_results, axis=0)
        isolate_summaries = pd.concat(isolate_summaries, axis=0)

        genes_affected = self.genes_affected(isolate_phenotypes)

        pointfinder_results_summary = self.pointfinder_summary(pointfinder_results)

        combinbined_presence_absence = pd.merge(
            genes_affected,
            pointfinder_results_summary,
            on="isolate_id",
        )

        return (
            pointfinder_results,
            resfinder_results,
            isolate_summaries,
            combinbined_presence_absence,
        )

    def collect(self):
        self.logger.info(f"Collecting results from {len(self.isolate_dirs)} isolates..")

        if self.output_dir is None:
            self.output_dir = self.RESFINDER_dir
        else:
            if os.path.isdir(self.output_dir) is False:
                os.mkdir(self.output_dir)

        (
            pointfinder_results,
            resfinder_results,
            isolate_summaries,
            combinbined_presence_absence,
        ) = self.collect_all_results()

        if pointfinder_results is None:
            self.logger.info("No results found.")
            return

        pointfinder_results.to_csv(
            os.path.join(self.output_dir, "pointfinder_results.tsv"),
            sep="\t",
            index=False,
        )
        resfinder_results.to_csv(
            os.path.join(self.output_dir, "resfinder_results.tsv"),
            sep="\t",
            index=False,
        )
        isolate_summaries.to_csv(
            os.path.join(self.output_dir, "isolate_summaries.tsv"),
            sep="\t",
            index=False,
        )

        combinbined_presence_absence.to_csv(
            os.path.join(self.output_dir, "combined_presence_absence.tsv"),
            sep="\t",
            index=False,
        )

        self.logger.info(f"Results written to {self.output_dir}")
