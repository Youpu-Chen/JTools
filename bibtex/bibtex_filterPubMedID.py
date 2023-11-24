#!/usr/bin/env python
# coding: utf-8
import bibtexparser
from bibtexparser.bwriter import BibTexWriter
import sys


def read_pubmed_ids(file_path):
    """Reads a file containing PubMed IDs (one per line) and returns them as a list."""
    with open(file_path, 'r') as file:
        return [line.strip() for line in file if line.strip()]


def retrieve_records_by_pubmed_id(bibtex_file, pubmed_ids):
    with open(bibtex_file) as bibtex_file:
        bib_database = bibtexparser.load(bibtex_file)

    targeted_records = []
    for entry in bib_database.entries:
        if 'pmid' in entry and entry['pmid'] in pubmed_ids:
            targeted_records.append(entry)

    return targeted_records


def write_bibtex_file(records, output_file):
    """Writes the given BibTeX records to a new file."""
    db = bibtexparser.bibdatabase.BibDatabase()
    db.entries = records

    writer = BibTexWriter()
    with open(output_file, 'w') as bibtex_file:
        bibtex_file.write(writer.write(db))


def main():
    pubmed_ids = read_pubmed_ids(sys.argv[1])  # Replace with your targeted PubMed ID list
    bibtex_file = sys.argv[2]  # Replace with your BibTeX file path
    records = retrieve_records_by_pubmed_id(bibtex_file, pubmed_ids)
    write_bibtex_file(records, sys.argv[3])



