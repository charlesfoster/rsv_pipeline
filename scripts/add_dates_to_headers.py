#!/usr/bin/env python3
"""
Append collection dates to FASTA headers.

Given:
  * a FASTA file whose headers start with a sequence ID
    (optionally followed by ‘|something’), and
  * a CSV file with two columns:  seq_id,date

the script rewrites each header as

    >seq_id|date

if a matching `seq_id` – exact string match with the part of the header
_before the first “|”_ – is found in the CSV.  Headers without a match are
left unchanged.

Usage
-----
    add_dates.py -i sequences.fasta -m dates.csv -o sequences_dated.fasta
"""
import csv
import argparse
from pathlib import Path


def load_dates(csv_path: Path) -> dict:
    """Return {seq_id: date} mapping from the two-column CSV."""
    mapping = {}
    with csv_path.open() as fh:
        rdr = csv.reader(fh)
        header = next(rdr, None)
        # tolerate header row – assume the first cell contains "seq_id"
        if header and header[0].lower().startswith("seq"):
            pass  # header already consumed
        else:     # data row already read – put it back
            if header:
                mapping[header[0]] = header[1]
        for seq_id, date in rdr:
            mapping[seq_id.strip()] = date.strip()
    return mapping


def rewrite_fasta(in_fa: Path, out_fa: Path, date_map: dict) -> None:
    with in_fa.open() as fin, out_fa.open("w") as fout:
        for line in fin:
            if line.startswith(">"):
                hdr = line[1:].rstrip()           # drop leading '>'
                seq_id = hdr.split("|", 1)[0]     # part before first '|'
                date   = date_map.get(seq_id)
                if date:
                    # remove any existing date after the first “|”
                    new_hdr = f"{seq_id}|{date}"
                    fout.write(f">{new_hdr}\n")
                else:
                    fout.write(line)              # unchanged
            else:
                fout.write(line)


def main() -> None:
    ap = argparse.ArgumentParser(description="Append dates to FASTA headers.")
    ap.add_argument("-i", "--input_fasta", required=True,
                    type=Path, help="Input FASTA")
    ap.add_argument("-m", "--metadata_csv", required=True,
                    type=Path, help="CSV with seq_id,date columns")
    ap.add_argument("-o", "--output_fasta", required=True,
                    type=Path, help="Output FASTA with dates in headers")
    args = ap.parse_args()

    date_map = load_dates(args.metadata_csv)
    rewrite_fasta(args.input_fasta, args.output_fasta, date_map)


if __name__ == "__main__":
    main()