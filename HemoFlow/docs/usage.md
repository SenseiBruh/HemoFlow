# Usage

Commands below run from the repository folder. On Windows the launchers change
to their own directory before starting Python.

## Viewer options

```powershell
.\START_HemoFlow.bat --defaults --cells 80 --fps 45
.\START_HemoFlow.bat --interactive --fps 45
.\START_HemoFlow.bat --defaults --mode single --cells 80 --bpm 72 --particles-off
```

`--interactive` prompts for physical inputs before the viewer starts. A JSON
configuration can also be supplied with `--config path\to\settings.json`.
`settings.example.json` documents the available fields. Command-line overrides
take precedence over loaded values. `--help` lists supported options.

The desktop grid-size ceiling limits interactive configurations. The N160/L150
headless exception is implemented only by the grid study wrapper.

## Recording

The headless interface accepts `--seconds`, `--cycles`, or `--continuous` with
`--record`. `--record-every` is a simulated-time sampling interval in seconds.
The viewer's record controls attach the same scalar recorder to the worker.

The recording directory contains time-series CSV, metadata JSON, and geometry
epochs. Headless final output includes numerical fields, settings, metrics, and
a geometry report. Headless fluid runs do not advect particles; particle
response quantities in those exports are reference parameters.

```powershell
.\.venv\Scripts\python.exe analyze_recording.py "exports\demo\raw\*_timeseries.csv" --csv-output "exports\demo\cycle_summary.csv"
```

The pattern must identify exactly one recording. Cycle summaries use
time-weighted interpolation and explicitly mark partial cycles. Add
`--complete-only` to exclude partial cycles; a run shorter than one cardiac
period will not produce a complete-cycle summary. CSV files open in Excel;
they can then be saved as an Excel workbook.

## Research studies

Study commands are documented in [studies/README.md](../studies/README.md).
They use their own source snapshots. Long-run results should be generated with
the same package that defines the experiment and its analysis criteria.
