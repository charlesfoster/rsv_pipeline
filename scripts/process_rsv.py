#!/usr/bin/env python3
# usage: python process_rsv.py -i in.fasta -d dates.csv -o out.fasta -c meta.csv
import re, csv, argparse
from pathlib import Path
import country_converter as coco

# ------------------ REGEX PATTERNS -------------------------------------------
GISAID_RE = re.compile(
    r'^(?P<virus>[^|]+)\|(?P<epi>EPI_ISL_[^|]+)\|(?P<date>\d{4}(?:-\d{1,2}(?:-\d{1,2})?)?)$'
)
DATE_RE   = re.compile(r'\d{4}(?:-\d{1,2}(?:-\d{1,2})?)?$')
RSV_PLAIN = re.compile(r'^RSV-(\d+)_')

# ------------------ HELPERS ---------------------------------------------------
def normalize_date(s):
    p = s.split('-')
    if len(p) == 3:
        return f"{p[0]}-{p[1].zfill(2)}-{p[2].zfill(2)}"
    if len(p) == 2:
        return f"{p[0]}-{p[1].zfill(2)}"
    return p[0]

def load_dates(csv_path):
    m = {}
    with open(csv_path) as fh:
        rdr = csv.reader(fh)
        first = next(rdr, None)
        if first and not first[0].lower().startswith('seq'):
            fh.seek(0)                       # no header
            rdr = csv.reader(fh)
        else:
            pass                             # header skipped
        for sid, date in rdr:
            m[sid.strip()] = normalize_date(date.strip())
    return m

def fix_header(h):
    if '|EPI_ISL_' in h or h.startswith('2024'):
        return h
    m = RSV_PLAIN.match(h)
    if m:
        return f"Peds2022RSV-{m.group(1)}"
    return h.split(' ', 1)[0]

def parse_header(h_line):
    full = h_line[1:].strip()                # drop '>'
    m = GISAID_RE.match(full)
    if m:                                    # ---------- GISAID ----------
        virus, epi, date = m.group('virus'), m.group('epi'), normalize_date(m.group('date'))
        parts   = virus.split('/')
        subtype = f"RSV{parts[1]}" if len(parts) > 1 else ""
        country = parts[2].replace('_', ' ') if len(parts) > 2 else ""
        return dict(
            seq_id            = f"{subtype}|{epi}|{date}",
            full_seq_header   = full,
            subtype           = subtype,
            virus_name        = virus,
            gisaid_accession  = epi,
            date_of_collection= date,
            country           = country,
            continent         = coco.convert(names=country, to="continent"),
            un_region         = coco.convert(names=country, to="UNregion"),
        )
    # -------------------- non-GISAID -----------------------------------------
    date = ""
    left, *right = full.rsplit('|', 1)
    if right and DATE_RE.fullmatch(right[0].strip()):
        date = normalize_date(right[0].strip())
    return dict(
        seq_id            = full,
        full_seq_header   = full,
        subtype           = "",
        virus_name        = "",
        gisaid_accession  = "",
        date_of_collection= date,
        country           = "",
        continent         = "",
        un_region         = "",
    )

# ------------------ MAIN ------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input",      required=True)
    ap.add_argument("-d", "--dates_csv",  required=True)
    ap.add_argument("-o", "--out_fasta",  required=True)
    ap.add_argument("-c", "--out_csv",    required=True)
    args = ap.parse_args()

    date_map = load_dates(Path(args.dates_csv))

    records = []
    header  = None
    seq_buf = []

    with open(args.input) as fin, open(args.out_fasta, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                if header is not None:                       # flush previous
                    rec = parse_header(header)
                    rec["sequence"] = "".join(seq_buf)
                    records.append(rec)
                    fout.write(f">{rec['seq_id']}\n{rec['sequence']}\n")
                raw   = line[1:].rstrip()
                fixed = fix_header(raw)
                sid   = fixed.split("|", 1)[0]
                if sid in date_map:
                    fixed = f"{sid}|{date_map[sid]}"
                header  = f">{fixed}"
                seq_buf = []
            else:
                seq_buf.append(line.strip())

        if header is not None:                               # final record
            rec = parse_header(header)
            rec["sequence"] = "".join(seq_buf)
            records.append(rec)
            fout.write(f">{rec['seq_id']}\n{rec['sequence']}\n")

    fields = [
        "seq_id", "full_seq_header", "subtype", "virus_name",
        "gisaid_accession", "date_of_collection", "country",
        "continent", "un_region"
    ]
    with open(args.out_csv, "w", newline="") as csvf:
        w = csv.DictWriter(csvf, fieldnames=fields)
        w.writeheader()
        for r in records:
            w.writerow({k: r[k] for k in fields})

if __name__ == "__main__":
    main()