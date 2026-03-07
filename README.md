# VFDistiller — Variant Fusion Distiller

A bioinformatics desktop tool for processing, converting, and annotating genetic variant data. Supports VCF, gVCF, 23andMe raw data, and FASTA — without pysam/bcftools/samtools (Windows-compatible).

![Variant Fusion - Main View](README/screenshots/main_view.png)

## Features

- **Multi-Format Import** — VCF, gVCF, 23andMe (.txt), FASTA (.fa/.fasta)
- **Automatic Build Detection** — GRCh37 / GRCh38 from header, contigs, or RSID positions
- **Multi-Source Annotation** — gnomAD, MyVariant.info, Ensembl VEP, ALFA, TOPMed, AlphaGenome
- **INFO Recycling** — Existing VCF annotations are reused
- **Filtering** — AF threshold, CADD score, Variant Impact, ClinSig, gene lists, FILTER=PASS, Read Depth
- **Export** — CSV, Excel, PDF, annotated VCF (filtered or complete)
- **GUI** — ttkbootstrap interface with System Tray, progress indicator, themes
- **Performance** — Optional Cython hot-path (5x overall speedup), SQLite batch writes, async HTTP via aiohttp
- **Background Maintenance** — Automatic re-fetching of missing annotations during idle
- **Multilingual** — German and English (JSON-based translations)

## Prerequisites

- Python 3.10+
- Windows 10/11 (primarily tested), Linux/macOS experimental

### Installation

```bash
# Install dependencies
pip install -r requirements.txt

# Optional: Cython acceleration (requires C compiler)
pip install cython
cd cython_hotpath
python setup.py build_ext --inplace
cd ..
```

### Genome References (optional, for FASTA validation)

The genome references (GRCh37/GRCh38) must be downloaded separately (~3 GB per build):

```bash
# GRCh37
wget https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz

# GRCh38
wget https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

Place the files in the project directory. On first launch, a `.fai` index is automatically generated.

### gnomAD LightDB (optional)

For fast offline AF lookups, the gnomAD LightDB can be downloaded. The tool offers a download dialog on first launch. Alternatively:

```bash
python "Get gnomAD DB light.py"
```

## Usage

### Launch GUI

```bash
python Variant_Fusion_pro_V17.py
```

Or on Windows:

```
START.bat
```

### Workflow

1. **Open file** — Select VCF, gVCF, 23andMe text file, or FASTA
2. **Check build** — Automatically detected, can be manually overridden
3. **Pipeline runs** — Variants are parsed, annotated, and filtered
4. **Results** — Table view with sortable columns, double-click opens external databases
5. **Export** — Export as CSV, Excel, PDF, or annotated VCF

### Configuration

On first launch, `variant_fusion_settings.json` is created from the template `variant_fusion_settings.json.example`. Key settings:

| Setting | Description | Default |
|---|---|---|
| `af_threshold` | Allele frequency threshold | 0.007 |
| `include_none` | Show variants without AF | false |
| `cadd_highlight_threshold` | CADD score highlighting | 22.0 |
| `stale_days` | Days until AF refresh | 200 |
| `alphagenom_key` | Google AlphaGenome API key | (empty) |
| `quality_settings` | VCF record-level filter | see example |

### API Keys

- **AlphaGenome**: Requires a Google AI API key. Enter in `variant_fusion_settings.json` under `alphagenom_key` and `api_settings.phase6_ag.alphagenom.api_key`.
- **NCBI**: Optional for higher rate limits. Enter under `api_settings.global.ncbi_api_key`.

## Dependencies

### Core (required)

| Package | License | Purpose |
|---|---|---|
| requests | Apache 2.0 | HTTP requests |
| psutil | BSD | CPU/Memory monitoring |
| Pillow | PIL License | Icon/Image processing |
| intervaltree | Apache 2.0 | Genomic intervals |
| ttkbootstrap | MIT | Modern GUI themes |
| pystray | MIT | System Tray icon |
| aiohttp | Apache 2.0 | Async HTTP fetching |
| scipy | BSD | Statistics |

### Optional

| Package | License | Purpose |
|---|---|---|
| openpyxl | MIT | Excel export |
| reportlab | BSD | PDF export |
| numpy | BSD | Array operations |
| biopython | Biopython License | Sequence alignment |
| pyfaidx | MIT | FASTA indexing |
| cython | Apache 2.0 | Hot-path compilation |

## Cython Acceleration

Optional C-compiled hot-paths for critical operations:

| Module | Speedup | Function |
|---|---|---|
| `vcf_parser.pyx` | 8x | VCF line parsing |
| `af_validator.pyx` | 100x | AF validation |
| `key_normalizer.pyx` | 25x | Variant key normalization |
| `fasta_lookup.pyx` | 100x | FASTA sequence lookup |

Overall pipeline speedup: ~5x (50k variants: 15 min -> 3 min).

If Cython is not installed, Python fallbacks are used automatically.

## Project Structure

```
VFDistiller/
├── Variant_Fusion_pro_V17.py .... Main program (GUI + Pipeline)
├── requirements.txt ............. Python dependencies
├── variant_fusion_settings.json.example . Configuration template
├── VFDistiller.spec ............. PyInstaller build configuration
├── START.bat .................... Windows quick-start
│
├── cython_hotpath/ .............. Optional Cython modules
│   ├── __init__.py .............. CythonAccelerator main class
│   ├── vcf_parser.pyx .......... VCF parsing
│   ├── af_validator.pyx ......... AF validation
│   ├── key_normalizer.pyx ....... Key normalization
│   ├── fasta_lookup.pyx ......... FASTA lookup
│   ├── setup.py ................. Build script
│   └── test_performance.py ...... Benchmarks
│
├── data/annotations/ ............ Gene annotation data
│   ├── GRCh37.gtf.gz ........... Ensembl gene annotations
│   └── GRCh38.gtf.gz
│
├── locales/
│   └── translations.json ........ Translations (de/en)
│
├── ICO/ICO.ico .................. App icon
│
├── lightdb_index_worker.py ...... gnomAD LightDB background indexing
├── translator.py ................ Translation engine
├── translator_patch.py .......... Translation patches
├── manage_translations.py ....... Translation management
├── Get gnomAD DB light.py ....... gnomAD download tool
├── test_performance.py .......... Performance tests
│
├── ARCHITECTURE.md .............. Developer documentation
└── README/ ...................... Extended documentation & licenses
    └── licenses/
        ├── LICENSE.txt .......... Main license (English)
        ├── LICENSE.de.txt ....... Main license (German)
        └── THIRD_PARTY_LICENSES.txt . Third-party licenses
```

## License

**VFDistiller License v1.0** — Free to use, modification allowed, no resale. See [LICENSE](LICENSE) for details.

- Use for research, education, and personal purposes: **allowed**
- Modification and adaptation of source code: **allowed**
- Redistribution within your own organization: **allowed**
- Resale or commercial redistribution: **prohibited**
- This license applies to V17.x — successor versions may have different terms

The software is not medically validated and must not be used for clinical diagnoses or therapeutic decisions.

Third-party libraries are subject to their respective licenses (MIT, BSD, Apache 2.0). See `README/licenses/THIRD_PARTY_LICENSES.txt`.

> **Windows Store:** A pre-packaged version with additional features (Cython acceleration, offline database) will be available in the Microsoft Store soon.

## Version

V17.0 — Current production version (March 2026).

---

🇩🇪 [Deutsche Version](README.de.md)
