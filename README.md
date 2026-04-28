# BTO Data Pipeline

Data processing pipeline for the BTO instrument, transforming raw spacecraft telemetry into calibrated science-ready data products across four processing levels.

---

## Table of Contents

- [Pipeline Overview](#pipeline-overview)
- [Data Levels](#data-levels)
- [Getting Started](#getting-started)
- [Installation](#installation)
- [Usage](#usage)
- [Repository Structure](#repository-structure)
- [Contributing](#contributing)
- [Open Decisions](#open-decisions)

---

## Pipeline Overview

Raw telemetry received from the spacecraft is ingested, verified, and progressively calibrated through the following stages:

```
Raw Telemetry → L0 (Binary) → L1a (Raw FITS) → L1b (Calibrated FITS) → L2 (Science FITS)
                                                                      ↓
                                                              influxDB + Grafana
```

All pipeline code lives under version control on this server and is publicly available. Calibration files are hosted on HEASARC via CALDB.

---

## Data Levels

### L0 — Binary (Raw Layer)
Raw binary packages received from the spacecraft, sorted by date (`yyyy/mm/dd`).

| Package Type | Description |
|---|---|
| D6 | Spectrum / Histograms |
| D7 | Event packages |
| D8 | Housekeeping |
| D9 | Registers |
| F6 | Bootloader HK |

**Processing steps:**
- Package verification (checksum)
- Sort files by date

---

### L1a — FITS (Raw Layer)
Parsed and split FITS files, converted from binary with GPS→UTC time conversion and mission elapsed time computed.

**Output files:**

| File Pattern | Contents | Cadence |
|---|---|---|
| `csYYMMDDbto.hk` | Housekeeping + Bootloader | Hourly |
| `csYYMMDDbto.lc` | Histograms | Hourly |
| *(light curves)* | Light curves | Hourly |
| `csYYMMDDXXXt_bto.evt` | Events / photon list (TTE) | Per event (`XXX` = fraction of day, `t` = trigger ID) |

---

### L1b — FITS (Calibrated Layer)
Calibrated FITS files with ADC values converted to physical units. Hosted on HEASARC; calibration files served from CALDB.

**Calibrations applied:**

| Data Type | Calibration |
|---|---|
| Housekeeping | ADC → physical units |
| Histograms | ADC → energy + deadtime |
| Light curves | Energy bands |
| Photon list | ADC → energy + deadtime |

**Auxiliary inputs required:**
- BTO HK calibration *(link TBD)*
- BTO Energy calibration *(link TBD)*
- BTO Deadtime calibration
- COSI position information + orientation *(TBD)*

---

### L2 — FITS (Science Layer)
Reconstructed good photon data. Data products must conform to Gamma-ray Data Tools standards.

**Auxiliary inputs:**
- Response matrix from lab setup

---

### influxDB + Grafana
Time-stamped raw and calibrated data (housekeeping + histograms) are ingested into an influxDB database. A Grafana dashboard provides real-time monitoring *(public availability TBD)*.

---

## Getting Started

> **Placeholder** — detailed quickstart instructions will be added once the environment is stable.

You will need:
- Python ≥ 3.10
- `astropy`, `numpy`, `fitsio` *(full list in `requirements.txt`)*
- Access to the HEASARC CALDB for calibration files

```bash
# Clone the repository
git clone https://github.com/<org>/bto-pipeline.git
cd bto-pipeline

# Install dependencies
pip install -r requirements.txt
```

---

## Usage

> **Placeholder** — command-line interface and example runs will be documented here.

```bash
# Example: Process a day of raw telemetry through L0 → L1a
python run_pipeline.py --input /data/raw/2024/01/15 --level l1a

# Example: Apply calibration to produce L1b products
python run_pipeline.py --input /data/l1a/2024/01/15 --level l1b
```

---

## Repository Structure

> **Placeholder** — update once the directory layout is finalised.

```
bto-pipeline/
├── l0/                  # L0 binary processing
├── l1a/                 # L1a FITS conversion
├── l1b/                 # L1b calibration
├── l2/                  # L2 science products
├── caldb/               # Local calibration file handling
├── utils/               # Shared utilities (time conversion, checksum, etc.)
├── tests/               # Unit and integration tests
├── docs/                # Additional documentation
└── README.md
```

---

## Contributing

### Branch Strategy

This repository uses a **feature branch workflow**:

- `main` is the stable, reviewed branch. Direct pushes are not permitted.
- All new work is developed on a dedicated branch and merged via a Pull Request (PR).
- Only the repository maintainer may merge PRs into `main`, following code review.

### Branch Naming

Use the following prefixes to keep the branch list readable:

| Prefix | Use for |
|---|---|
| `feature/` | New functionality (e.g. `feature/l1b-energy-calibration`) |
| `fix/` | Bug fixes (e.g. `fix/gps-utc-conversion`) |
| `docs/` | Documentation only (e.g. `docs/update-readme`) |
| `refactor/` | Code cleanup with no functional change |
| `test/` | Adding or improving tests |

### Workflow Step-by-Step

1. **Create a branch** from the latest `main`:
   ```bash
   git checkout main && git pull
   git checkout -b feature/your-feature-name
   ```

2. **Develop and commit** your changes with clear, descriptive commit messages.

3. **Push your branch** and open a Pull Request on GitHub:
   ```bash
   git push origin feature/your-feature-name
   ```

4. **Fill in the PR template** — describe what changed, reference any relevant issues, and confirm local testing.

5. **Request review** from the maintainer. Address any feedback with follow-up commits on the same branch.

6. **Merge** — once approved, the maintainer merges the PR into `main`. Prefer *Squash and merge* for a clean history on small features, or *Merge commit* to preserve detailed history on larger efforts.

7. **Delete your branch** after merging to keep the repo tidy.

### Pull Request Checklist

Before requesting review, confirm:

- [ ] Code runs without errors locally
- [ ] New functionality includes tests (or tests are noted as a follow-up)
- [ ] Docstrings/comments updated where relevant
- [ ] No hardcoded paths or credentials
- [ ] FITS output validated against expected format
- [ ] Update the version number of the pipeline 

### Code Review Guidelines

Reviewers should check for:
- Correctness of calibration logic
- FITS header compliance
- Time system handling (GPS vs UTC vs MET)
- Consistent file naming conventions (see Data Levels above)

---

## Open Decisions

The following items are currently unresolved and tracked here until closed:

| Level | Question |
|---|---|
| L0 | Do we append packages? |
| L0 | Do we sort by time? (BTO data may not arrive time-sorted) |
| L1a | FITS Header format — port Carolyn's email |
| L1a | Mirror BTO ADC housekeeping into both data streams? |
| L1a | Should there be a filler for time periods with no data? |
| L1a | Do we archive register readbacks? |
| L1b | Config file format — FITS extensions, FITS header, or standalone file? |
| L1b | Visualization of data approach |
| L1b | BTO HK calibration link (missing) |
| L1b | BTO Energy calibration link (missing) |
| L1b | COSI position + orientation as auxiliary input? |
| L2 | Confirm conformance requirements with Gamma-ray Data Tools |
| influxDB | Document table schema and ingestion details |
| Grafana | Should the dashboard be publicly available? |

---

## License

This project is licensed under the **Apache License, Version 2.0**.

You are free to use, modify, and distribute this software for any purpose, including commercial use, provided that you include the required attribution notice. Contributors are also granted an explicit patent license under Apache 2.0.

See the [LICENSE](./LICENSE) file for the full license text, or visit [apache.org/licenses/LICENSE-2.0](https://www.apache.org/licenses/LICENSE-2.0).


