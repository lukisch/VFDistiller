# VFDistiller Migration Concept v1
## Zerlegung des Monolithen in Module

**Stand:** 16.01.2026
**Ausgangslage:** Variant_Fusion_pro_V13.py (~20.000+ Zeilen, "God-Script")
**Ziel:** Wartbare, testbare Paket-Struktur

---

## 1. Ziel-Architektur: Layered Architecture

Trennung in Schichten:
- **Daten** (DB)
- **Logik** (Pipeline/Distiller)
- **Netzwerk** (API)
- **Darstellung** (GUI)

---

## 2. Neue Ordnerstruktur

```
variant_fusion/
│
├── main.py                    # Einstiegspunkt (Startet App)
├── config.py                  # Zentrale Konfiguration (Config-Klasse)
├── constants.py               # Konstanten (Queries, Regex, Pfade)
│
├── core/                      # Kern-Logik & Manager
│   ├── __init__.py
│   ├── flags.py               # Flag_and_Options_Manager
│   ├── quality.py             # QualityManager
│   ├── filters.py             # MainFilterGate, CodingFilter, AfNoneTreatment
│   └── progress.py            # PipelineProgress
│
├── database/                  # Datenhaltung
│   ├── __init__.py
│   ├── db_manager.py          # VariantDB Klasse
│   ├── buffer.py              # VCFBuffer, EmitQueue
│   └── lightdb.py             # LightDBGnomADManager
│
├── pipeline/                  # Die Verarbeitungs-Engine
│   ├── __init__.py
│   ├── distiller.py           # Distiller (Orchestrator)
│   ├── scanner.py             # VCF Scanning Logik (_phase_vcf_scan)
│   ├── fetcher_pipeline.py    # Streaming Phasen (AF, Full-Anno)
│   └── background.py          # BackgroundMaintainer, BackofficeCrawler
│
├── network/                   # Externe Kommunikation
│   ├── __init__.py
│   ├── async_fetcher.py       # MyVariant async Logik
│   ├── controller.py          # AFFetchController
│   └── resilience.py          # CircuitBreaker, ThroughputTuner
│
├── bio/                       # Bioinformatische Helper & Parser
│   ├── __init__.py
│   ├── fasta.py               # FastaValidator
│   ├── genes.py               # GeneAnnotator
│   ├── vcf_parser.py          # parse_vcf_records_smart, etc.
│   ├── converters.py          # 23andMe, FastQ Konverter
│   └── alphagenome.py         # AlphaGenomeScorer
│
├── ui/                        # Grafische Oberfläche
│   ├── __init__.py
│   ├── app.py                 # Hauptklasse App (Tkinter root)
│   ├── windows.py             # QualitySettingsDialog, etc.
│   ├── components.py          # Widgets (ScrolledFrame, etc.)
│   └── export.py              # Export-Logik (VCF, Excel, PDF)
│
└── utils/                     # Allgemeine Hilfsfunktionen
    ├── __init__.py
    ├── logger.py              # MultiSinkLogger
    ├── formatting.py          # fmt_eta, safe_float, genotype helpers
    └── system.py              # Threading helpers, mmap utils
```

---

## 3. Detaillierte Klassen-Zuordnung

### core/ (Das Gehirn der Einstellungen)

| Klasse | Zieldatei | Hinweis |
|--------|-----------|---------|
| Flag_and_Options_Manager | core/flags.py | Muss isoliert sein (keine GUI-Importe) |
| QualityManager | core/quality.py | |
| MainFilterGate | core/filters.py | |
| CodingFilter | core/filters.py | |
| AfNoneTreatmentManager | core/filters.py | |
| FetchStatusManager | core/filters.py | |
| PipelineProgress | core/progress.py | |

### database/ (Das Gedächtnis)

| Klasse | Zieldatei | Hinweis |
|--------|-----------|---------|
| VariantDB | database/db_manager.py | SQL-Statements auslagern |
| VCFBuffer | database/buffer.py | |
| EmitQueue | database/buffer.py | Technisch ein Buffer zur GUI |
| LightDBGnomADManager | database/lightdb.py | |

### bio/ (Die Werkzeuge)

| Klasse/Funktion | Zieldatei | Hinweis |
|-----------------|-----------|---------|
| FastaValidator | bio/fasta.py | |
| ShiftDetectionStats | bio/fasta.py | |
| GeneAnnotator | bio/genes.py | |
| AlphaGenomeScorer | bio/alphagenome.py | |
| convert_23andme_to_vcf | bio/converters.py | |
| FASTQmap | bio/converters.py | |
| parse_vcf_records | bio/vcf_parser.py | |
| parse_vcf_records_smart | bio/vcf_parser.py | |
| parse_vcf_records_mmap | bio/vcf_parser.py | |
| is_gzipped | bio/vcf_parser.py | |

### pipeline/ (Der Prozess)

| Klasse | Zieldatei | Hinweis |
|--------|-----------|---------|
| Distiller | pipeline/distiller.py | Orchestrator, delegiert an Scanner/Fetcher |
| _phase_vcf_scan Logik | pipeline/scanner.py | Eigene Klasse VcfScanner |
| _phase_af_fetch_streaming | pipeline/fetcher_pipeline.py | |
| _phase_full_annotation_streaming | pipeline/fetcher_pipeline.py | |
| BackgroundMaintainer | pipeline/background.py | |
| BackofficeCrawler | pipeline/background.py | |
| VCFMigrationsdienst | pipeline/background.py | |

### network/ (Die Kommunikation)

| Klasse/Funktion | Zieldatei | Hinweis |
|-----------------|-----------|---------|
| AFFetchController | network/controller.py | |
| ThroughputTuner | network/resilience.py | |
| CircuitBreaker | network/resilience.py | |
| mv_fetch | network/async_fetcher.py | |
| gnomad_fetch_async | network/async_fetcher.py | |

### ui/ (Das Gesicht)

| Klasse | Zieldatei | Hinweis |
|--------|-----------|---------|
| App | ui/app.py | Export-Methoden auslagern! |
| QualitySettingsDialog | ui/windows.py | |
| Export-Funktionen | ui/export.py | Aus App auslagern |

### utils/ (Hilfsfunktionen)

| Klasse/Funktion | Zieldatei |
|-----------------|-----------|
| MultiSinkLogger | utils/logger.py |
| fmt_eta | utils/formatting.py |
| safe_float | utils/formatting.py |
| genotype helpers | utils/formatting.py |

---

## 4. Migrationsplan (Phasenweise)

### Phase 1: Basis und Utils (Geringstes Risiko)

1. Erstelle `utils/logger.py`, `utils/formatting.py`, `config.py`
2. Verschiebe entsprechende Klassen/Funktionen
3. Passe Importe im Hauptskript an
4. **Test:** Anwendung muss noch laufen

### Phase 2: Bio-Module und Datenbank

1. Extrahiere `VariantDB` nach `database/` (kritisch!)
2. Extrahiere `FastaValidator` und `GeneAnnotator` nach `bio/`
3. Diese sind relativ isoliert
4. **Test:** VCF-Import muss funktionieren

### Phase 3: Core Manager (Lösen der Abhängigkeiten)

1. Extrahiere `Flag_and_Options_Manager` und `QualityManager`
2. **WICHTIG:** Diese Module dürfen NICHT `App` importieren!
3. Die `App` importiert sie und setzt die Werte
4. Das löst viele zirkuläre Abhängigkeiten

### Phase 4: Die Pipeline zerlegen (Der schwierige Teil)

1. Löse `Distiller` aus der `App`
2. Statt `self.app` an Distiller zu übergeben: Callback-Funktionen nutzen
3. Trenne `_phase_vcf_scan` in eigene Klasse `VcfScanner`
4. **Test:** Kompletter Pipeline-Durchlauf

### Phase 5: GUI Cleanup

1. `App` Klasse enthält nur noch GUI-Code
2. Alle Logik delegiert an entsprechende Module
3. Export-Funktionen in `ui/export.py`

---

## 5. Abhängigkeits-Hierarchie (Import-Reihenfolge)

```
Level 0 (Keine Projekt-Abhängigkeiten):
├── config.py
├── constants.py
└── utils/*.py

Level 1 (Nur Level 0):
├── database/db_manager.py
├── database/lightdb.py
└── bio/*.py

Level 2 (Level 0 + 1):
├── core/*.py
├── network/*.py
└── database/buffer.py

Level 3 (Level 0 + 1 + 2):
└── pipeline/*.py

Level 4 (Alles):
└── ui/*.py

Level 5 (Entry Point):
└── main.py
```

**Regel:** Ein Modul darf nur Module aus niedrigeren Levels importieren!

---

## 6. Zirkuläre Import-Vermeidung

### Das Problem

```python
# pipeline.py
from ui import App  # Distiller braucht App

# ui.py
from pipeline import Distiller  # App braucht Distiller

# → ImportError!
```

### Die Lösung: Callback-Pattern

**Alt (im Distiller):**
```python
self.app.progress.set_phase("vcf_scan")
```

**Neu (im Distiller):**
```python
# pipeline/distiller.py
class Distiller:
    def __init__(self, db, progress_tracker, ...):
        self.progress = progress_tracker  # Instanz von PipelineProgress
        # Keine Referenz auf "app" nötig!

    def run(self):
        self.progress.set_phase("vcf_scan")
```

**In ui/app.py:**
```python
from core.progress import PipelineProgress
from pipeline.distiller import Distiller

class App(ttk.Window):
    def __init__(self):
        self.progress_tracker = PipelineProgress()
        self.distiller = Distiller(
            db=self.db,
            progress_tracker=self.progress_tracker
        )

    def _refresh_ui(self):
        # UI holt sich Daten vom Tracker, nicht umgekehrt
        status = self.progress_tracker.get_status()
        self.progressbar.set(status['percent'])
```

---

## 7. Service Router (Optional, für schrittweise Migration)

Für eine risikoarme Migration kann ein Router-Pattern verwendet werden:

```python
# core/router.py
class ServiceRouter:
    """Zentraler Router für Legacy/Neu-Umschaltung"""

    def __init__(self):
        # Feature Flags
        self.use_new_vcf_parser = False
        self.use_new_db_writer = False

    def parse_vcf(self, path):
        if self.use_new_vcf_parser:
            from bio.vcf_parser import VcfScanner
            return VcfScanner().parse(path)
        else:
            # Legacy
            from Variant_Fusion_pro_V13 import parse_vcf_records_smart
            return parse_vcf_records_smart(path)

# Globale Instanz
router = ServiceRouter()
```

**Vorteile:**
- Risikominimierung (Fallback auf alten Code)
- Schrittweise Migration
- A/B Testing möglich

---

## 8. Einfacher Ansatz: Dateien trennen ohne Aufruf-Änderungen

### Reihenfolge für sicheres Trennen

**Schritt 1: Die "Unabhängigen" (Level 0)**
```python
# config.py - Config-Klasse
# utils.py - MultiSinkLogger, fmt_eta, safe_float
# bio_tools.py - FastaValidator, GeneAnnotator
```
Diese importieren nichts aus dem Projekt.

**Schritt 2: Die Datenbank (Level 1)**
```python
# database.py - VariantDB, VCFBuffer
```
Braucht nur Config und Logger.

**Schritt 3: Die Logik (Level 2)**
```python
# engine.py - Distiller, MainFilterGate, QualityManager
```
Braucht DB (Level 1) und Tools (Level 0).

**Schritt 4: Die GUI (Level 3)**
```python
# gui.py - App, QualitySettingsDialog
```
Darf alles importieren.

### Resultierende main.py

```python
# main.py

# Level 0
from config import Config
from utils import MultiSinkLogger

# Level 1
from database import VariantDB

# Level 2
from engine import Distiller, Flag_and_Options_Manager

# Level 3
from gui import App

if __name__ == "__main__":
    logger = MultiSinkLogger(...)
    db = VariantDB(...)
    app = App()
    app.mainloop()
```

---

## 9. Kritische Abhängigkeit zum Kappen

**Problem:** `VariantDB`/`VCFBuffer` greift auf `self.distiller.emit_queue` zu.

**Lösung:**
- Entferne harte Kopplung
- Nutze Callback-Funktion oder Event-System
- DB sollte "dumm" sein (nur Daten speichern)

**Suche im Code nach:**
```python
self.distiller  # innerhalb VariantDB oder VCFBuffer
self.app        # innerhalb Distiller oder Phasen
```

Diese Referenzen durch Callbacks ersetzen!

---

## 10. Roadmap für Router-Migration

| Phase | Funktion | Router-Methode | Warum zuerst? |
|-------|----------|----------------|---------------|
| 1 | VCF Parsing | `router.scan_vcf(path)` | Größter Bottleneck, isolierte Logik |
| 2 | Datenbank | `router.save_variants(batch)` | Kritisch für Stabilität |
| 3 | API Fetching | `router.fetch_annotations(keys)` | Kapselt aiohttp/requests Komplexität |
| 4 | Export | `router.export_data(format)` | Gut isolierbar |

---

## 11. Zusammenfassung

**Empfohlener Ansatz:**

1. **Kurzfristig:** Dateien trennen nach Level-Hierarchie
2. **Mittelfristig:** Router-Pattern für riskante Teile (Parser, DB-Writer)
3. **Langfristig:** Vollständige Modul-Struktur mit Callback-Entkopplung

**Kritischer Erfolgsfaktor:**
Die Abhängigkeit `DB → Distiller` und `Distiller → App` kappen durch:
- Callbacks statt direkter Referenzen
- Event-System für Updates
- Dependency Injection

**Testbarkeit als Ziel:**
Pipeline soll ohne GUI-Fenster testbar sein!
