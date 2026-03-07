# Store Listing — VFDistiller

## Deutsch

### Kurzbeschreibung (max 100 Zeichen)
Genetische Varianten analysieren, annotieren und filtern — VCF, gVCF, 23andMe

### Beschreibung (max 10.000 Zeichen)
VFDistiller (Variant Fusion Distiller) ist ein bioinformatisches Desktop-Tool zur Verarbeitung, Konvertierung und Annotation genetischer Variantendaten.

FEATURES:
- Multi-Format-Import: VCF, gVCF, 23andMe-Rohdaten (.txt), FASTA (.fa/.fasta)
- Automatische Build-Erkennung: GRCh37 und GRCh38 aus Header, Contigs oder RSID-Positionen
- Multi-Source-Annotation: gnomAD, MyVariant.info, Ensembl VEP, ALFA, TOPMed, AlphaGenome
- INFO-Recycling: Vorhandene VCF-Annotationen werden intelligent wiederverwendet
- Umfangreiche Filterung: AF-Schwelle, CADD-Score, Variant Impact, ClinSig, Genlisten, FILTER=PASS, Read Depth
- Vielseitiger Export: CSV, Excel, PDF, annotiertes VCF (gefiltert oder vollstaendig)
- Moderne GUI: ttkbootstrap-Oberflaeche mit System-Tray, Fortschrittsanzeige und Themes
- Performance: Optionaler Cython-Hotpath (5x Gesamt-Speedup), SQLite-Batch-Writes, async HTTP
- Hintergrund-Wartung: Automatisches Nachladen fehlender Annotationen im Leerlauf
- Mehrsprachig: Deutsch und Englisch

FUER WEN:
- Bioinformatiker und Genetiker, die Variantendaten effizient analysieren moechten
- Forschende, die VCF-Dateien aus verschiedenen Quellen zusammenfuehren und annotieren wollen
- 23andMe/AncestryDNA-Nutzer, die ihre Rohdaten wissenschaftlich auswerten moechten
- Alle, die eine Windows-native Alternative zu pysam/bcftools/samtools suchen

VORTEILE:
- Keine pysam/bcftools/samtools noetig — laeuft nativ unter Windows
- Offline-faehig dank gnomAD LightDB fuer schnelle Allele-Frequency-Lookups
- Automatischer Download fehlender Genomreferenzen und Datenbanken
- Transparente Ergebnisse: Tabellenansicht mit sortierbaren Spalten, Doppelklick oeffnet externe Datenbanken

KOSTENLOSE BASIS-VERSION
Diese Version ist kostenlos und enthaelt alle Kern-Features.

### Schluesselwoerter
VCF, Varianten, Genetik, Bioinformatik, gnomAD, Annotation, 23andMe, FASTA, Genom, Analyse, Filter, CADD, ClinVar

### Kategorie
Photo & Video

---

## English

### Short Description (max 100 chars)
Analyze, annotate and filter genetic variants — VCF, gVCF, 23andMe support

### Description (max 10,000 chars)
VFDistiller (Variant Fusion Distiller) is a bioinformatics desktop tool for processing, converting and annotating genetic variant data.

FEATURES:
- Multi-Format Import: VCF, gVCF, 23andMe raw data (.txt), FASTA (.fa/.fasta)
- Automatic Build Detection: GRCh37 and GRCh38 from header, contigs or RSID positions
- Multi-Source Annotation: gnomAD, MyVariant.info, Ensembl VEP, ALFA, TOPMed, AlphaGenome
- INFO Recycling: Existing VCF annotations are intelligently reused
- Comprehensive Filtering: AF threshold, CADD score, variant impact, ClinSig, gene lists, FILTER=PASS, read depth
- Versatile Export: CSV, Excel, PDF, annotated VCF (filtered or complete)
- Modern GUI: ttkbootstrap interface with system tray, progress indicators and themes
- Performance: Optional Cython hot-path (5x overall speedup), SQLite batch writes, async HTTP
- Background Maintenance: Automatic re-fetching of missing annotations during idle time
- Multilingual: German and English

FOR WHOM:
- Bioinformaticians and geneticists who want to efficiently analyze variant data
- Researchers who need to merge and annotate VCF files from multiple sources
- 23andMe/AncestryDNA users who want to scientifically evaluate their raw data
- Anyone looking for a Windows-native alternative to pysam/bcftools/samtools

BENEFITS:
- No pysam/bcftools/samtools required — runs natively on Windows
- Offline-capable with gnomAD LightDB for fast allele frequency lookups
- Automatic download of missing genome references and databases
- Transparent results: sortable table view, double-click opens external databases

FREE BASE VERSION
This version is free and includes all core features.

### Keywords
VCF, variants, genetics, bioinformatics, gnomAD, annotation, 23andMe, FASTA, genome, analysis, filter, CADD, ClinVar

### Category
Photo & Video
