#!/usr/bin/env python3

"""
Stage 4 Taxonomic Annotation Engine

Resolves BLAST accession identifiers from the Stage 3 species composition
table into NCBI taxonomy information and produces an enriched CSV output.
"""


from __future__ import annotations

import json
import logging
import re
import sys
import time
import xml.etree.ElementTree as ET

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional


import pandas as pd
import requests

from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


###############################################################################
# STEP 1 Configuration
###############################################################################

VERSION = "1.0.0"

NCBI_EMAIL = "rudie.kauhanen@sussex.ac.uk"

NCBI_TOOL = "sussex_edna_pipeline"

REQUEST_TIMEOUT = 30

REQUEST_DELAY = 0.34

USER_AGENT = (
    f"{NCBI_TOOL}/{VERSION} "
    f"({NCBI_EMAIL})"
)

###############################################################################
# Paths
###############################################################################

BASE_DIR = Path(__file__).resolve().parent

RUNS_DIR = BASE_DIR / "runs"


if not RUNS_DIR.exists():
    raise FileNotFoundError(
        f"Missing runs directory:\n{RUNS_DIR}"
    )


RUNS = sorted(
    path for path in RUNS_DIR.iterdir()
    if path.is_dir()
)


if not RUNS:
    raise RuntimeError(
        "No pipeline runs found."
    )


CURRENT_RUN = RUNS[-1]


SPECIES_DIRS = sorted(
    path for path in CURRENT_RUN.iterdir()
    if (
        path.is_dir()
        and path.name.startswith("final_species_list_")
    )
)


if not SPECIES_DIRS:
    raise RuntimeError(
        "No Stage 3 species directory found."
    )


STAGE3_DIR = SPECIES_DIRS[-1]


INPUT_CSV = (
    STAGE3_DIR /
    "master_species_composition.csv"
)


OUTPUT_CSV = (
    STAGE3_DIR /
    "master_species_composition_taxonomy.csv"
)


CACHE_FILE = (
    STAGE3_DIR /
    "taxonomy_cache.json"
)


LOG_FILE = (
    STAGE3_DIR /
    "taxonomy.log"
)


if not INPUT_CSV.exists():
    raise FileNotFoundError(
        f"Missing input CSV:\n{INPUT_CSV}"
    )

ACCESSION_PATTERN = re.compile(
    r"([A-Z]{1,4}_?\d+\.\d+)"
)

###############################################################################
# Logging
###############################################################################

LOGGER = logging.getLogger("taxonomy_engine")

LOGGER.setLevel(logging.INFO)

LOGGER.handlers.clear()


handler = logging.StreamHandler(sys.stdout)

handler.setFormatter(
    logging.Formatter(
        "%(asctime)s | %(levelname)s | %(message)s"
    )
)

LOGGER.addHandler(handler)


file_handler = logging.FileHandler(
    LOG_FILE,
    encoding="utf-8"
)

file_handler.setFormatter(
    logging.Formatter(
        "%(asctime)s | %(levelname)s | %(message)s"
    )
)

LOGGER.addHandler(file_handler)

###############################################################################
# STEP 2 Taxonomy Data Model
###############################################################################

@dataclass(slots=True)
class TaxonomyRecord:

    accession: str

    taxid: str = ""

    scientific_name: str = ""

    rank: str = ""

    kingdom: str = ""

    phylum: str = ""

    class_name: str = ""

    order: str = ""

    family: str = ""

    genus: str = ""

    species: str = ""

    status: str = "Unresolved"

###############################################################################
# Cache Manager
###############################################################################

class CacheManager:
    """
    Persistent accession -> taxonomy cache.
    """

    def __init__(self, path: Path):

        self.path = path

        self.records: dict[str, TaxonomyRecord] = {}

        self.load()

    def load(self):

        if not self.path.exists():

            LOGGER.info(
                "No taxonomy cache found."
            )

            return

        with open(
            self.path,
            "r",
            encoding="utf-8"
        ) as file:

            data = json.load(file)

        for accession, record in data.items():

            self.records[accession] = TaxonomyRecord(
                **record
            )

        LOGGER.info(
            "Loaded %d cached records.",
            len(self.records)
        )

    def save(self):

        with open(
            self.path,
            "w",
            encoding="utf-8"
        ) as file:

            json.dump(

                {
                    accession: asdict(record)

                    for accession, record
                    in sorted(self.records.items())

                },

                file,

                indent=2

            )

    def get(
        self,
        accession: str
    ) -> Optional[TaxonomyRecord]:

        return self.records.get(
            accession
        )

    def put(
        self,
        record: TaxonomyRecord
    ):
        self.records[
            record.accession
        ] = record

###############################################################################
# STEP 3 NCBI Client
###############################################################################

class NCBIClient:

    BASE_URL = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    )

    def __init__(self):

        self.session = requests.Session()

        retry = Retry(
            total=5,
            backoff_factor=1,
            status_forcelist=[429, 500, 502, 503, 504],
            allowed_methods=["GET"]
        )

        adapter = HTTPAdapter(
            max_retries=retry
        )

        self.session.mount(
            "https://",
            adapter
        )

        self.session.headers.update(
            {
                "User-Agent": USER_AGENT
            }
        )


    def request(
        self,
        endpoint: str,
        **params
    ):

        params.update(
            {
                "tool": NCBI_TOOL,
                "email": NCBI_EMAIL
            }
        )

        response = self.session.get(
            f"{self.BASE_URL}/{endpoint}",
            params=params,
            timeout=REQUEST_TIMEOUT
        )

        response.raise_for_status()

        time.sleep(
            REQUEST_DELAY
        )

        return response



    def accession_to_taxid(
        self,
        accessions: list[str]
    ) -> dict[str, str]:

        response = self.request(
            "esummary.fcgi",
            db="nuccore",
            id=",".join(accessions),
            retmode="json"
        )

        data = response.json()

        result = data.get(
            "result",
            {}
        )

        mapping = {}

        for uid in result.get(
            "uids",
            []
        ):

            record = result.get(
                uid,
                {}
            )

            accession = record.get(
                "accessionversion"
            )

            taxid = record.get(
                "taxid"
            )

            if accession and taxid:

                mapping[accession] = str(taxid)

        return mapping



    def get_taxonomy(
        self,
        taxids: list[str]
    ) -> dict[str, TaxonomyRecord]:

        response = self.request(
            "efetch.fcgi",
            db="taxonomy",
            id=",".join(taxids),
            retmode="xml"
        )

        root = ET.fromstring(
            response.text
        )

        records = {}

        for taxon in root.findall(
            ".//Taxon"
        ):

            taxid = taxon.findtext(
                "TaxId",
                ""
            )

            name = taxon.findtext(
                "ScientificName",
                ""
            )

            rank = taxon.findtext(
                "Rank",
                ""
            )

            lineage = {
                "kingdom": "",
                "phylum": "",
                "class": "",
                "order": "",
                "family": "",
                "genus": ""
            }

            lineage_xml = taxon.find(
                "LineageEx"
            )

            if lineage_xml is not None:

                for ancestor in lineage_xml.findall(
                    "Taxon"
                ):

                    ancestor_rank = ancestor.findtext(
                        "Rank",
                        ""
                    ).lower()

                    ancestor_name = ancestor.findtext(
                        "ScientificName",
                        ""
                    )

                    if ancestor_rank == "superkingdom":

                        lineage["kingdom"] = ancestor_name

                    elif ancestor_rank in lineage:

                        lineage[ancestor_rank] = ancestor_name

            records[taxid] = TaxonomyRecord(
                accession="",
                taxid=taxid,
                scientific_name=name,
                rank=rank,
                kingdom=lineage["kingdom"],
                phylum=lineage["phylum"],
                class_name=lineage["class"],
                order=lineage["order"],
                family=lineage["family"],
                genus=lineage["genus"],
                species=name if rank.lower() == "species" else "",
                status="Resolved"
            )

        return records



    def resolve(
        self,
        accessions: list[str]
    ) -> dict[str, TaxonomyRecord]:

        accession_taxids = self.accession_to_taxid(
            accessions
        )

        taxonomy = {}

        if accession_taxids:
            taxonomy = self.get_taxonomy(
        list(accession_taxids.values())
    )

        resolved = {}

        for accession, taxid in accession_taxids.items():

            record = taxonomy.get(
                taxid
            )

            if record:

                record.accession = accession

                resolved[accession] = record

        return resolved
    
###############################################################################
# STEP 4 Taxonomy Resolver
###############################################################################

class TaxonomyResolver:
    def __init__(self, ncbi_client: NCBIClient, cache_manager: CacheManager):
        self.ncbi = ncbi_client
        self.cache = cache_manager

    @staticmethod
    def normalise_accession(accession: str) -> str:
        if not accession:
            return ""

        accession = str(accession).strip()

        match = ACCESSION_PATTERN.search(accession)

        return match.group(1) if match else accession

    def resolve(self, accessions: list[str]) -> dict[str, TaxonomyRecord]:
        cleaned = {
            self.normalise_accession(acc)
            for acc in accessions
            if acc
        }

        LOGGER.info(
            "Unique accession identifiers: %d",
            len(cleaned)
        )

        resolved = {}
        missing = []

        for accession in cleaned:
            record = self.cache.get(accession)

            if record:
                resolved[accession] = record
            else:
                missing.append(accession)

        LOGGER.info(
            "Cache hits: %d",
            len(resolved)
        )

        LOGGER.info(
            "NCBI queries required: %d",
            len(missing)
        )

        if missing:
            ncbi_results = self.ncbi.resolve(missing)

            for accession, record in ncbi_results.items():
                resolved[accession] = record
                self.cache.put(record)

        LOGGER.info(
            "Total taxonomy records resolved: %d",
            len(resolved)
        )

        return resolved
    
###############################################################################
# STEP 5 Taxonomy Report Builder
###############################################################################

class TaxonomyReportBuilder:

    def __init__(self, resolver: TaxonomyResolver):
        self.resolver = resolver

    def build_report(self, input_csv: Path, output_csv: Path):

        LOGGER.info("Loading %s", input_csv)

        df = pd.read_csv(input_csv)

        accession_column = next(
            (
                col for col in [
                    "Accession",
                    "accession",
                    "AccessionID"
                ]
                if col in df.columns
            ),
            None
        )

        if accession_column is None:
            raise ValueError(
                "No accession column found."
            )

        accessions = (
            df[accession_column]
            .dropna()
            .astype(str)
            .tolist()
        )

        taxonomy = self.resolver.resolve(accessions)

        taxonomy_rows = []

        for accession in df[accession_column]:

            clean = TaxonomyResolver.normalise_accession(
                accession
            )

            record = taxonomy.get(clean)

            if record:

                taxonomy_rows.append(
                    {
                        "NCBI_TaxID": record.taxid,
                        "NCBI_Scientific_Name": record.scientific_name,
                        "NCBI_Rank": record.rank,
                        "Kingdom": record.kingdom,
                        "Phylum": record.phylum,
                        "Class": record.class_name,
                        "Order": record.order,
                        "Family": record.family,
                        "Genus": record.genus,
                        "Species": record.species,
                        "Taxonomy_Status": "Resolved"
                    }
                )

            else:

                taxonomy_rows.append(
                    {
                        "NCBI_TaxID": "",
                        "NCBI_Scientific_Name": "",
                        "NCBI_Rank": "",
                        "Kingdom": "",
                        "Phylum": "",
                        "Class": "",
                        "Order": "",
                        "Family": "",
                        "Genus": "",
                        "Species": "",
                        "Taxonomy_Status": "Unresolved"
                    }
                )

        result = pd.concat(
            [
                df.reset_index(drop=True),
                pd.DataFrame(taxonomy_rows)
            ],
            axis=1
        )

        result.to_csv(
            output_csv,
            index=False
        )

        LOGGER.info(
            "Written taxonomy report: %s",
            output_csv
        )

        return result
    
###############################################################################
# STEP 6 Main Execution
###############################################################################

def main():

    LOGGER.info(
        "Starting taxonomy annotation."
    )

    cache = CacheManager(
        CACHE_FILE
    )

    ncbi = NCBIClient()

    resolver = TaxonomyResolver(
        ncbi,
        cache
    )

    builder = TaxonomyReportBuilder(
        resolver
    )

    builder.build_report(
        INPUT_CSV,
        OUTPUT_CSV
    )

    cache.save()

    LOGGER.info(
        "Taxonomy annotation complete."
    )


if __name__ == "__main__":
    main()