# Research studies

These packages preserve the numerical implementations used for the recent
domain and grid studies. Each has a local `source_frozen` directory and checksum
manifest. Those intentional copies support reproducibility and are independent
of the desktop application's root modules.

| Package | Comparison | Solver |
| --- | --- | --- |
| [length_confirmation](length_confirmation/README.md) | N80 at 120 and 150 mm | Frozen section-impedance outlet |
| [grid_refinement](grid_refinement/README.md) | N80, N120, conditional N160 at 150 mm | Frozen section-impedance outlet |

The numerical code, plans, tolerances, and historical setup-check files are
preserved from the supplied packages. The README files have been rewritten.
New runs will therefore record a different package-documentation hash.
Existing studies must be resumed with their original package and source receipt.

From the repository folder, inspect either plan without running the solver:

```powershell
.\.venv\Scripts\python.exe studies\length_confirmation\run_length_confirmation.py --plan-only
.\.venv\Scripts\python.exe studies\grid_refinement\run_grid.py --plan-only
```

Study-specific unit tests run in separate processes to keep their frozen module
imports independent of the desktop application:

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s studies\length_confirmation -p "test_*.py"
.\.venv\Scripts\python.exe -m unittest discover -s studies\grid_refinement -p "test_*.py"
```

The grid runner retains its conditional N160 behavior and has no wall-time
limit. Fine-grid runs can require multiple days. `--plan-only` inspects the
sequence without starting a batch; runtime estimates are machine-specific.

Older study packages remain in the accompanying original source archive. Their
results require their original configuration and analysis definitions.
