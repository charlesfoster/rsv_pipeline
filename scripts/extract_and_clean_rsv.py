#!/usr/bin/env python3
import re
import csv
import argparse
import country_converter as coco

# ---- Regex for GISAID-style headers -----------------------------------------
GISAID_RE = re.compile(
    r"""^
        (?P<virus>[^|]+)       \|
        (?P<epi>EPI_ISL_[^|]+) \|
        (?P<date>\d{4}(?:-\d{1,2}(?:-\d{1,2})?)?)
    $""",
    re.VERBOSE,
)

# ---- Loose date regex -------------------------------------------------------
DATE_RE = re.compile(r"\d{4}(?:-\d{1,2}(?:-\d{1,2})?)?$")

# ---- Utilities --------------------------------------------------------------
def normalize_date(date_str: str) -> str:
    parts = date_str.split("-")
    if len(parts) == 3:
        y, m, d = parts
        return f"{y}-{m.zfill(2)}-{d.zfill(2)}"
    if len(parts) == 2:
        y, m = parts
        return f"{y}-{m.zfill(2)}"
    return parts[0]

# ---- Header parser ----------------------------------------------------------
def parse_header(header_line: str) -> dict:
    full = header_line[1:].strip()
    m = GISAID_RE.match(full)

    if m:                                    # ---------- GISAID ----------
        virus = m.group("virus")
        epi   = m.group("epi")
        date  = normalize_date(m.group("date"))

        parts = virus.split("/")
        subtype = f"RSV{parts[1]}" if len(parts) > 1 else ""
        country = parts[2].replace("_", " ") if len(parts) > 2 else ""
        continent = coco.convert(names=country, to="continent")
        un_region = coco.convert(names=country, to="UNregion")
        seq_id = f"{subtype}|{epi}|{date}"

        return {
            "seq_id": seq_id,
            "full_seq_header": full,
            "subtype": subtype,
            "virus_name": virus,
            "gisaid_accession": epi,
            "date_of_collection": date,
            "country": country,
            "continent": continent,
            "un_region": un_region,
        }

    # ------------------------- non-GISAID ------------------------------------
    date = ""
    left, *right = full.rsplit("|", 1)          # look ONLY at text after last '|'
    if right and DATE_RE.fullmatch(right[0].strip()):
        date = normalize_date(right[0].strip())

    return {
        "seq_id": full,
        "full_seq_header": full,
        "subtype": "",
        "virus_name": "",
        "gisaid_accession": "",
        "date_of_collection": date,
        "country": "",
        "continent": "",
        "un_region": "",
    }

# ---- Main -------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Extract RSV metadata and rewrite FASTA headers."
    )
    ap.add_argument("-i", "--input",     required=True, help="Input FASTA")
    ap.add_argument("-c", "--out_csv",   required=True, help="Output metadata CSV")
    ap.add_argument("-o", "--out_fasta", required=True, help="Output cleaned FASTA")
    args = ap.parse_args()

    records = []
    with open(args.input) as fin, open(args.out_fasta, "w") as fout:
        header = None
        seq_buf = []
        for line in fin:
            if line.startswith(">"):
                if header is not None:
                    rec = parse_header(header)
                    rec["sequence"] = "".join(seq_buf)
                    records.append(rec)
                    fout.write(f">{rec['seq_id']}\n{rec['sequence']}\n")
                header, seq_buf = line, []
            else:
                seq_buf.append(line.strip())

        if header is not None:
            rec = parse_header(header)
            rec["sequence"] = "".join(seq_buf)
            records.append(rec)
            fout.write(f">{rec['seq_id']}\n{rec['sequence']}\n")

    # ---- CSV ----------------------------------------------------------------
    fields = [
        "seq_id",
        "full_seq_header",
        "subtype",
        "virus_name",
        "gisaid_accession",
        "date_of_collection",
        "country",
        "continent",
        "un_region",
    ]
    with open(args.out_csv, "w", newline="") as csvf:
        writer = csv.DictWriter(csvf, fieldnames=fields)
        writer.writeheader()
        for rec in records:
            writer.writerow({k: rec[k] for k in fields})

if __name__ == "__main__":
    main()
    