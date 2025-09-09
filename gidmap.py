#!/usr/bin/env python3
"""
gidmap.py — A tiny, fast CLI tool to convert gene identifiers

Backed by MyGene.info (https://mygene.info) via the `mygene` Python package.
Works for many species (default: human, 9606) and many ID types.

Features
--------
- Convert from one ID type (e.g., symbol) to one or many target types (e.g., entrezgene, ensembl.gene, uniprot).
- Batch converts from a file, STDIN, or command-line list.
- TSV/CSV/JSON output; prints to STDOUT by default (or use -o to write a file).
- Returns best hit by default; can expand to all hits with --all.
- Handles unmapped IDs gracefully and preserves input order.
- Optional local HTTP cache for speed/retries.

Installation (recommended: pipx)
--------------------------------
1) Ensure pipx is installed

   pip install pipx

2) In your cloned repo (with pyproject.toml), run

   pipx install .

3) Use it from anywhere as `gidmap`.

Minimal usage
-------------
echo -e "TP53\nBRCA1\nEGFR" | gidmap --from-type symbol --to entrezgene ensembl.gene --species human

Example output (TSV by default):
input_id\tfrom_type\tentrezgene\tensembl.gene\tmatched_symbol\tmatch_type
TP53\tsymbol\t7157\tENSG00000141510\tTP53\texact
BRCA1\tsymbol\t672\tENSG00000012048\tBRCA1\texact
EGFR\tsymbol\t1956\tENSG00000146648\tEGFR\texact

Supported identifier types
--------------------------
FROM (use with --from-type / --scopes): symbol, entrezgene, ensembl.gene, ensembl.transcript, ensembl.protein, uniprot, refseq, alias, mgd, wormbase, flybase, zfin, tair, etc.

TO (use with --to / --fields): entrezgene, symbol, name, ensembl.gene, ensembl.transcript, ensembl.protein, uniprot, refseq.rna, refseq.protein, taxonomy, genomic_pos, map_location, etc.

Notes
-----
- Species can be a common name ("human", "mouse") or NCBI taxid (9606, 10090). Defaults to 9606.
- MyGene scopes/fields are documented here: https://docs.mygene.info/en/latest/doc/query_service.html
- For fully offline mapping, you would need a local mapping table; this script uses the live API for freshness.
"""
from __future__ import annotations

import sys
import argparse
import csv
import json
from typing import Iterable, List, Dict, Any, Optional

try:
    import mygene  # type: ignore
except Exception as e:
    sys.stderr.write(
        "[gidmap] ERROR: The 'mygene' package is required. Install with: pip install mygene\n"
    )
    raise

# Optional caching for reliability/speed
try:
    import requests_cache  # type: ignore
except Exception:
    requests_cache = None  # type: ignore


def read_ids(args: argparse.Namespace) -> List[str]:
    ids: List[str] = []
    if args.ids:
        ids = list(
            dict.fromkeys([s.strip() for s in args.ids if s.strip()])
        )  # unique preserve order
    elif args.input and args.input != "-":
        with open(args.input, "r", encoding=args.encoding) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if args.column:
                    # If a delimited file with a column index (1-based)
                    parts = line.split(args.delimiter)
                    idx = args.column - 1
                    if 0 <= idx < len(parts):
                        ids.append(parts[idx].strip())
                else:
                    ids.append(line)
    else:
        # STDIN
        for line in sys.stdin:
            line = line.strip()
            if not line:
                continue
            if args.column:
                parts = line.split(args.delimiter)
                idx = args.column - 1
                if 0 <= idx < len(parts):
                    ids.append(parts[idx].strip())
            else:
                ids.append(line)
    return ids


def normalize_species(spec: str) -> str | int:
    # Accept simple aliases
    s = spec.strip().lower()
    aliases = {
        "human": 9606,
        "homo sapiens": 9606,
        "mouse": 10090,
        "mus musculus": 10090,
        "rat": 10116,
        "rattus norvegicus": 10116,
        "zebrafish": 7955,
        "danio rerio": 7955,
        "fly": 7227,
        "drosophila melanogaster": 7227,
        "yeast": 559292,
        "saccharomyces cerevisiae": 559292,
        "ecoli": 511145,
        "escherichia coli": 511145,
    }
    if s in aliases:
        return aliases[s]
    # If it's an int, assume taxid
    try:
        return int(s)
    except ValueError:
        return spec  # mygene also accepts species names


def flatten_field(doc: Dict[str, Any], key: str) -> Any:
    """Safely pull nested fields with dots in `key` (e.g., 'ensembl.gene').
    If list, join with ';'. If dict, attempt common keys.
    """
    cur: Any = doc
    for part in key.split("."):
        if isinstance(cur, list):
            # flatten each list item dict by accessing the field if present
            cur = [item.get(part) if isinstance(item, dict) else item for item in cur]
        elif isinstance(cur, dict):
            cur = cur.get(part)
        else:
            cur = None
        if cur is None:
            break
    # At this point cur might be list/dict/str/int
    if isinstance(cur, list):
        # Deduplicate while preserving order
        seen = set()
        vals = []
        for x in cur:
            if x is None:
                continue
            if isinstance(x, dict):
                # pick an obvious ID-like value if present
                for k in ("gene", "transcript", "protein", "id", "accession"):
                    if k in x and x[k] is not None:
                        v = str(x[k])
                        if v not in seen:
                            seen.add(v)
                            vals.append(v)
                        break
            else:
                v = str(x)
                if v not in seen:
                    seen.add(v)
                    vals.append(v)
        return ";".join(vals) if vals else None
    if isinstance(cur, dict):
        # Common ID keys
        for k in ("gene", "transcript", "protein", "id", "accession"):
            if k in cur and cur[k] is not None:
                return cur[k]
        return json.dumps(cur, ensure_ascii=False)
    return cur


def convert_ids(
    ids: List[str],
    scopes: str,
    fields: List[str],
    species: str | int,
    allow_multiple: bool,
    chunksize: int = 1000,
    retries: int = 2,
    as_dict: bool = False,
) -> List[Dict[str, Any]]:
    mg = mygene.MyGeneInfo()

    results: List[Dict[str, Any]] = []

    def do_query(batch: List[str]) -> List[Dict[str, Any]]:
        return mg.querymany(
            batch,
            scopes=scopes,
            fields=fields + ["symbol", "preferred_symbol", "name"],
            species=species,
            returnall=False,
            as_dataframe=False,
            verbose=False,
            dotfield=True,
        )

    # MyGene returns results not strictly preserving order when multi-hits; we'll re-attach input_id
    for i in range(0, len(ids), chunksize):
        batch = ids[i : i + chunksize]
        attempt = 0
        while True:
            try:
                res = do_query(batch)
                break
            except Exception as e:
                attempt += 1
                if attempt > retries:
                    raise
        # Normalize each record
        by_input: Dict[str, List[Dict[str, Any]]] = {}
        for r in res:
            _in = r.get("query") or r.get("input")
            if _in is None:
                continue
            by_input.setdefault(str(_in), []).append(r)
        for q in batch:
            hits = by_input.get(q, [])
            if not hits:
                results.append(
                    {
                        "input_id": q,
                        "from_type": scopes,
                        "match_type": "not_found",
                    }
                )
                continue
            if allow_multiple:
                for h in hits:
                    results.append(_format_hit(q, scopes, fields, h))
            else:
                # pick highest score / best match
                hbest = max(hits, key=lambda d: d.get("_score", 0))
                results.append(_format_hit(q, scopes, fields, hbest))

    return results


def _format_hit(
    query: str, scopes: str, fields: List[str], hit: Dict[str, Any]
) -> Dict[str, Any]:
    row: Dict[str, Any] = {
        "input_id": query,
        "from_type": scopes,
    }
    # match type
    mt = "match"
    if hit.get("notfound"):
        mt = "not_found"
    elif hit.get("disambiguation"):
        mt = "ambiguous"
    elif hit.get("_score") is not None:
        mt = "exact" if hit.get("score") == hit.get("_score") else "best"
    row["match_type"] = mt

    # include a helpful symbol/name if available
    matched_symbol = hit.get("symbol") or hit.get("preferred_symbol")
    if matched_symbol:
        row["matched_symbol"] = matched_symbol
    if "name" in hit:
        row["name"] = hit["name"]

    # populate target fields
    for f in fields:
        row[f] = flatten_field(hit, f)
    return row


def write_output(
    rows: List[Dict[str, Any]], outfmt: str, outfile: Optional[str]
) -> None:
    if outfmt == "json":
        data = rows
        text = json.dumps(data, ensure_ascii=False, indent=2)
        if outfile:
            with open(outfile, "w", encoding="utf-8") as f:
                f.write(text + "\n")
        else:
            sys.stdout.write(text + "\n")
        return

    # For tabular outputs, construct fieldnames from union of keys preserving a sensible order
    base_cols = ["input_id", "from_type", "matched_symbol", "name", "match_type"]
    extra_cols = []
    seen = set(base_cols)
    for r in rows:
        for k in r.keys():
            if k not in seen and k not in base_cols:
                extra_cols.append(k)
                seen.add(k)
    fieldnames = base_cols + extra_cols

    dialect = "excel-tab" if outfmt == "tsv" else "excel"
    out = open(outfile, "w", encoding="utf-8", newline="") if outfile else sys.stdout
    writer = csv.DictWriter(out, fieldnames=fieldnames, dialect=dialect)
    writer.writeheader()
    for r in rows:
        writer.writerow({k: r.get(k, "") for k in fieldnames})
    if outfile:
        out.close()


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        prog="gidmap",
        description="Convert gene identifiers via MyGene.info",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "ids", nargs="*", help="IDs to convert (optional if --input or stdin used)"
    )
    p.add_argument(
        "--from-type",
        "--scopes",
        dest="scopes",
        required=True,
        help="Source identifier type (MyGene 'scopes'), e.g., symbol, entrezgene, ensembl.gene, uniprot, alias",
    )
    p.add_argument(
        "--to",
        "--fields",
        dest="fields",
        nargs="+",
        required=True,
        help="One or more target fields to retrieve, e.g., entrezgene symbol ensembl.gene uniprot refseq.rna",
    )
    p.add_argument("--species", default="9606", help="Species name or taxid")

    io = p.add_argument_group("Input options")
    io.add_argument(
        "-i", "--input", default="-", help="Input file path, or '-' for stdin"
    )
    io.add_argument(
        "--delimiter",
        default="\t",
        help="Delimiter for parsing a column from a file/stdin",
    )
    io.add_argument(
        "--column",
        type=int,
        help="1-based column index to read IDs from when parsing delimited text",
    )
    io.add_argument("--encoding", default="utf-8")

    outg = p.add_argument_group("Output options")
    outg.add_argument(
        "-o", "--output", help="Write results to this file (default: stdout)"
    )
    outg.add_argument("--format", choices=["tsv", "csv", "json"], default="tsv")
    outg.add_argument(
        "--all",
        action="store_true",
        help="Return all hits (may produce multiple rows per input)",
    )

    net = p.add_argument_group("Network/cache")
    net.add_argument(
        "--cache", metavar="PATH", help="Enable HTTP cache file (via requests-cache)"
    )

    args = p.parse_args(argv)

    # caching
    if args.cache and requests_cache is not None:
        requests_cache.install_cache(args.cache)

    species = normalize_species(args.species)

    # Read list of input IDs
    inputs = read_ids(args)
    if not inputs:
        sys.stderr.write(
            "[gidmap] No input IDs provided. Pass as arguments, --input file, or via stdin.\n"
        )
        return 2

    # Convert
    rows = convert_ids(
        inputs,
        scopes=args.scopes,
        fields=args.fields,
        species=species,
        allow_multiple=args.all,
    )

    # Write
    write_output(rows, outfmt=args.format, outfile=args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
