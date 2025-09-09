# gidmap

*A tiny, fast CLI to convert gene identifiers via [MyGene.info](https://mygene.info).*  
**Works across many species** (default: human, taxid 9606) and many ID types (NCBI/Ensembl/UniProt/RefSeq, aliases, and model-organism DBs).

---

## Key Features

* Convert from one ID type (e.g., `symbol`) to one or many target fields (e.g., `entrezgene`, `ensembl.gene`, `uniprot`).
* Batch convert from **file**, **STDIN**, or a **CLI list**.
* **TSV/CSV/JSON** output; prints to STDOUT by default (or `-o` to write a file).
* Returns **best hit** by default; optional `--all` to list all hits.
* Handles **unmapped IDs** cleanly and preserves **input order**.
* Optional local **HTTP cache** for speed/retries.

> Powered by the excellent [MyGene.info API](https://docs.mygene.info/).

---

## Installation (recommended: `pipx`)

```bash
# one-time setup (if you don't have pipx yet)
pip install pipx

# install gidmap from your cloned repo root
pipx install .

# now it's available on your PATH
gidmap --help
```

> Alternative installs
>
> * **venv**: `python -m venv .venv && source .venv/bin/activate && pip install -e .`
> * **conda/mamba**: create env and `pip install -e .` inside it
> * **Docker**: see below

---

## Quick Start

Convert human symbols to Entrez & Ensembl gene IDs:

```bash
echo -e "TP53\nBRCA1\nEGFR" | gidmap --from-type symbol --to entrezgene ensembl.gene --species human
```

**Example (TSV)**

```
input_id\tfrom_type\tentrezgene\tensembl.gene\tmatched_symbol\tmatch_type
TP53\tsymbol\t7157\tENSG00000141510\tTP53\texact
BRCA1\tsymbol\t672\tENSG00000012048\tBRCA1\texact
EGFR\tsymbol\t1956\tENSG00000146648\tEGFR\texact
```

Read IDs from a file (use column 2 in a TSV):

```bash
gidmap -i ids.tsv --column 2 --from-type symbol --to entrezgene symbol ensembl.gene --species 9606 -o out.tsv
```

Return **all** matching hits (not just best):

```bash
gidmap --from-type alias --to entrezgene symbol ensembl.gene --all -i aliases.txt
```

JSON output with a local HTTP cache:

```bash
gidmap --from-type ensembl.gene --to symbol entrezgene --species mouse \
       --format json --cache .gidmap_cache -i ens_ids.txt
```

---

## Supported Identifier Types

**Input scopes** (`--from-type`):

* **Symbols & aliases**: `symbol`, `alias`
* **NCBI**: `entrezgene`, `refseq` (`refseq.rna`, `refseq.protein`)
* **Ensembl**: `ensembl.gene`, `ensembl.transcript`, `ensembl.protein`
* **Protein**: `uniprot` (Swiss-Prot/TrEMBL)
* **Model organisms**: `mgd` (MGI), `wormbase`, `flybase`, `zfin`, `tair`, etc.

**Output fields** (`--to`):

* Basic: `entrezgene`, `symbol`, `name`, `alias`
* Ensembl: `ensembl.gene`, `ensembl.transcript`, `ensembl.protein`
* UniProt: `uniprot.Swiss-Prot`, `uniprot.TrEMBL`
* RefSeq: `refseq.rna`, `refseq.protein`
* Position/Mapping: `genomic_pos`, `map_location`
* Organism/Meta: `taxid`, `species`

> For the complete, up-to-date field catalog, check the MyGene.info docs:
> [https://docs.mygene.info/en/latest/doc/annotation\_service.html#available\_fields](https://docs.mygene.info/en/latest/doc/annotation_service.html#available_fields)

---

## CLI Reference

```
gidmap [IDS ...] --from-type <SCOPES> --to <FIELD1> [<FIELD2> ...] [options]

Options:
  --species <name|taxid>     Species name or NCBI taxid (default: 9606)
  -i, --input <path|->       Input file path, or '-' for stdin (default: '-')
      --delimiter <char>     Delimiter for parsing a column (default: '\t')
      --column <1-based>     Column index to read IDs from when parsing delimited text
      --encoding <codec>     Input encoding (default: utf-8)
  -o, --output <path>        Write results to this file (default: stdout)
      --format <tsv|csv|json> Output format (default: tsv)
      --all                  Return all hits (may produce multiple rows per input)
      --cache <path>         Enable HTTP cache (via requests-cache)
  -h, --help
```

**Species shortcuts**: you can pass common names (*human*, *mouse*, *rat*, *yeast*, *ecoli*, etc.) or numeric taxids (e.g., `9606`).

---

## Tips & Gotchas

* **Ambiguous symbols**: Some symbols are reused across species/gene families. Use a precise `--species` to improve matches. Use `--all` to inspect alternatives.
* **Multiple IDs per field**: Fields like `ensembl.transcript` or `uniprot` may return multiple values; `gidmap` joins them with `;` in tabular outputs.
* **Aliases**: Searching with `--from-type alias` is broad but can include historical/secondary names; confirm results with `matched_symbol`.
* **Reproducibility**: Use `--cache` for speed; for strict reproducibility, consider containerizing (see below).

---

## Docker (optional)

Create a `Dockerfile` in the repo root:

```dockerfile
FROM python:3.11-slim
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
RUN pip install --no-cache-dir mygene requests requests-cache
WORKDIR /app
COPY gidmap.py /app/gidmap.py
ENTRYPOINT ["python", "/app/gidmap.py"]
```

Build & run:

```bash
docker build -t gidmap:0.1 .
echo -e "TP53\nBRCA1" | docker run -i --rm gidmap:0.1 --from-type symbol --to entrezgene ensembl.gene --species human
```

## Acknowledgements

This tool relies on the **MyGene.info** service. Please cite their resources where appropriate and be mindful of API usage guidelines.

---

## License

MIT
