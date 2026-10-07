#!/usr/bin/env python3
"""
Validate and prepare finished consensus genomes for Nextclade and UShER.

One FASTA file per sample; the sample ID is the file name without extension
(e.g. SAMPLE001.fasta -> SAMPLE001). Original headers (e.g. SPAdes NODE_1_...) are
replaced by the sample ID. Every FASTA must match a samplesheet row and vice
versa - any problem stops the run before analysis.

Sequences are oriented to the reference with minimap2, since assemblies may be
reverse-complemented relative to NC_063383.

Outputs (in the working directory):
  {run}_combined_consensus.fasta   one record per sample, header = sample ID
  {run}_samplesheet_consensus.csv  samplesheet (';') with sample_id + tip_label
  {run}_input_check.tsv            per-sample input summary
"""
import argparse
import csv
import gzip
import re
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

FASTA_EXT = re.compile(r"\.(fasta|fa|fna)(\.gz)?$", re.IGNORECASE)
# Characters that break Newick tip labels or the VCF/SAM sample columns
UNSAFE_CHARS = re.compile(r"[\s:,();'\[\]]")
COMPLEMENT = str.maketrans("ACGTRYKMBDHVNacgtrykmbdhvn", "TGCAYRMKVHDBNtgcayrmkvhdbn")


def fail(problems):
    print("ERROR: consensus input validation failed:", file=sys.stderr)
    for p in problems:
        print(f"  - {p}", file=sys.stderr)
    sys.exit(1)


def read_fasta(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    records = []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                records.append([line[1:], []])
            elif records:
                records[-1][1].append(line)
    return [(header, "".join(seq).upper()) for header, seq in records]


def read_samplesheet(path):
    raw = Path(path).read_bytes()
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = raw.decode("latin-1")
    reader = csv.DictReader(text.splitlines(), delimiter=";")
    fieldnames = [f.strip() for f in reader.fieldnames]
    rows = []
    for row in reader:
        row = {k.strip(): (v or "").strip() for k, v in row.items() if k is not None}
        if any(row.values()):
            rows.append(row)
    return fieldnames, rows


def orientation(fasta_path, ref):
    """Return {seq_id: (plus_matches, minus_matches)} from a minimap2 PAF."""
    paf = subprocess.run(
        ["minimap2", "-c", "-x", "asm20", "--secondary=no", ref, fasta_path],
        check=True, capture_output=True, text=True,
    ).stdout
    matches = defaultdict(lambda: [0, 0])
    for line in paf.splitlines():
        cols = line.split("\t")
        matches[cols[0]][0 if cols[4] == "+" else 1] += int(cols[9])
    return matches


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--ref", required=True, help="reference FASTA used for orientation")
    ap.add_argument("--run-name", required=True)
    ap.add_argument("--fasta", nargs="+", required=True)
    args = ap.parse_args()
    run = args.run_name
    problems = []

    # --- FASTA files: ID from file name, exactly one record each ---------------
    samples = {}  # id -> dict
    ids_seen = Counter()
    for path in args.fasta:
        name = Path(path).name
        if not FASTA_EXT.search(name):
            problems.append(f"{name}: not a FASTA file name (.fasta/.fa/.fna, optionally .gz)")
            continue
        sid = FASTA_EXT.sub("", name)
        ids_seen[sid] += 1
        if UNSAFE_CHARS.search(sid):
            problems.append(f"{name}: sample ID '{sid}' contains whitespace or one of : , ( ) ; ' [ ]")
        records = read_fasta(path)
        if len(records) != 1:
            problems.append(
                f"{name}: contains {len(records)} sequences, expected exactly 1 "
                "(multi-contig assemblies must be scaffolded/merged first)"
            )
            continue
        header, seq = records[0]
        if not seq:
            problems.append(f"{name}: sequence is empty")
            continue
        samples[sid] = {"file": name, "header": header, "seq": seq}
    problems += [f"sample ID '{sid}' comes from {n} FASTA files" for sid, n in ids_seen.items() if n > 1]

    # --- Samplesheet ------------------------------------------------------------
    fieldnames, rows = read_samplesheet(args.samplesheet)
    id_col = "sample_id" if "sample_id" in fieldnames else "PrøveID" if "PrøveID" in fieldnames else None
    missing = [c for c in ("RunName", "SampleDate") if c not in fieldnames]
    if id_col is None:
        missing.insert(0, "PrøveID (or sample_id)")
    if missing:
        fail([f"samplesheet is missing required column(s): {', '.join(missing)}"])

    sheet_ids = Counter(r[id_col] for r in rows)
    problems += [f"samplesheet: sample ID '{sid}' appears {n} times" for sid, n in sheet_ids.items() if n > 1]
    problems += [f"samplesheet: row {i} has an empty {id_col}" for i, r in enumerate(rows, 2) if not r[id_col]]
    run_names = sorted({r["RunName"] for r in rows})
    if run_names != [run]:
        problems.append(f"samplesheet: RunName must be the single value '{run}', found {run_names}")
    for r in rows:
        try:
            datetime.strptime(r["SampleDate"], "%d.%m.%Y")
        except ValueError:
            problems.append(f"samplesheet: sample '{r[id_col]}' has invalid SampleDate '{r['SampleDate']}' (expected DD.MM.YYYY)")

    problems += [f"FASTA for '{sid}' ({samples[sid]['file']}) has no samplesheet row" for sid in samples if sid not in sheet_ids]
    problems += [f"samplesheet row '{sid}' has no FASTA file" for sid in sheet_ids if sid and sid not in ids_seen]

    if problems:
        fail(problems)

    # --- Orientation against the reference --------------------------------------
    tmp = "unoriented.fasta"
    with open(tmp, "w") as fh:
        for sid, s in samples.items():
            fh.write(f">{sid}\n{s['seq']}\n")
    matches = orientation(tmp, args.ref)
    unmapped = [sid for sid in samples if sum(matches.get(sid, (0, 0))) == 0]
    if unmapped:
        fail([f"'{sid}' does not align to the MPXV reference" for sid in unmapped])
    for sid, s in samples.items():
        plus, minus = matches[sid]
        s["revcomp"] = minus > plus
        if s["revcomp"]:
            s["seq"] = s["seq"].translate(COMPLEMENT)[::-1]

    # --- Outputs (samplesheet order) ----------------------------------------------
    order = [r[id_col] for r in rows]
    with open(f"{run}_combined_consensus.fasta", "w") as fh:
        for sid in order:
            seq = samples[sid]["seq"]
            fh.write(f">{sid}\n")
            fh.writelines(seq[i:i + 60] + "\n" for i in range(0, len(seq), 60))

    out_fields = ["sample_id" if f == id_col else f for f in fieldnames] + ["tip_label"]
    with open(f"{run}_samplesheet_consensus.csv", "w", newline="") as fh:
        writer = csv.writer(fh, delimiter=";")
        writer.writerow(out_fields)
        for r in rows:
            writer.writerow([r.get(f, "") for f in fieldnames] + [r[id_col]])

    with open(f"{run}_input_check.tsv", "w") as fh:
        fh.write("sample_id\tsource_file\toriginal_header\tlength\tpercent_N\treverse_complemented\n")
        for sid in order:
            s = samples[sid]
            pct_n = 100 * s["seq"].count("N") / len(s["seq"])
            fh.write(f"{sid}\t{s['file']}\t{s['header']}\t{len(s['seq'])}\t{pct_n:.2f}\t{'yes' if s['revcomp'] else 'no'}\n")

    flipped = sum(s["revcomp"] for s in samples.values())
    print(f"Prepared {len(samples)} sequences for run {run} ({flipped} reverse-complemented)")


if __name__ == "__main__":
    main()
