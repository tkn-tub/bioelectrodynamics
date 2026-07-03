# Running the ProDy frequency calculation on an HPC

This package converts the Colab notebook into a two-stage workflow:

1. **HPC**: compute ANM frequencies only, without plotting.
2. **Local PC**: copy the CSV results back and plot them locally.

The main calculation script is:

```bash
compute_prody_frequencies.py
```

It produces:

```text
results/frequencies_long.csv
results/summary.csv
results/frequencies.npz
results/run_manifest.json
```

The local plotting script is:

```bash
plot_frequencies_local.py
```

It reads `frequencies_long.csv` and writes PNG files.

---

## 1. Copy files to the HPC

From your local computer:

```bash
scp hpc_prody_package.zip YOUR_USERNAME@YOUR_HPC_ADDRESS:~/prody_job/
```

Then log in:

```bash
ssh YOUR_USERNAME@YOUR_HPC_ADDRESS
cd ~/prody_job
unzip hpc_prody_package.zip
```

---

## 2. Create the Python environment on the HPC

### Option A: venv + pip

Use this if your HPC has Python 3.12.

```bash
module avail python
module load Python/3.12

python -m venv ~/venvs/prody-hpc
source ~/venvs/prody-hpc/bin/activate

python -m pip install -U pip setuptools wheel
python -m pip install -r requirements-hpc.txt

python -c "import prody, numpy; print('ProDy', prody.__version__); print('NumPy', numpy.__version__)"
```

### Option B: conda or mamba

Use this if your HPC provides Anaconda/Miniconda/Mamba.

```bash
module avail conda
module load Miniconda3

conda env create -f environment.yml
conda activate prody-hpc

python -c "import prody, numpy; print('ProDy', prody.__version__); print('NumPy', numpy.__version__)"
```

If your HPC uses `mamba`, this is often faster:

```bash
mamba env create -f environment.yml
mamba activate prody-hpc
```

---

## 3. Download structure files

Many HPC compute nodes do not have internet access. Download the structures on the login node first:

```bash
source ~/venvs/prody-hpc/bin/activate
python compute_prody_frequencies.py --download-only --structures-dir structures --outdir results
```

For conda, activate your conda environment instead of the venv.

This creates a `structures/` directory containing the PDB files used by the computation.

---

## 4. Run the calculation

### Interactive test

Run a small test on the login node or an interactive compute session:

```bash
python compute_prody_frequencies.py \
  --structures-dir structures \
  --outdir results \
  --no-download \
  --n-modes 20
```

If that works, run all modes:

```bash
python compute_prody_frequencies.py \
  --structures-dir structures \
  --outdir results \
  --no-download \
  --n-modes all
```

### SLURM batch job

Edit `run_prody.slurm` and adjust the environment activation line. Then submit:

```bash
sbatch run_prody.slurm
```

Check the job:

```bash
squeue -u $USER
tail -f prody_JOBID.out
```

Replace `JOBID` with the number shown by `sbatch`.

---

## 5. Copy results back to your local PC

From your local computer:

```bash
scp -r YOUR_USERNAME@YOUR_HPC_ADDRESS:~/prody_job/results .
```

---

## 6. Plot locally

On your local PC:

```bash
python -m venv .venv-plot
```

Windows PowerShell:

```powershell
.\.venv-plot\Scripts\Activate.ps1
python -m pip install -U pip
python -m pip install -r requirements-local-plot.txt
python plot_frequencies_local.py --input results/frequencies_long.csv --outdir local_plots
```

macOS/Linux:

```bash
source .venv-plot/bin/activate
python -m pip install -U pip
python -m pip install -r requirements-local-plot.txt
python plot_frequencies_local.py --input results/frequencies_long.csv --outdir local_plots
```

The main plot will be:

```text
local_plots/vdos_overlay.png
```

Individual histograms will be in:

```text
local_plots/histograms/
```

---

## Notes

- The HPC calculation script intentionally does **not** import `matplotlib`, `seaborn`, `animation`, or `IPython.display`.
- The default frequency conversion is the one from your notebook: `frequency_THz = 0.225 * sqrt(eigenvalue)`.
- If the job runs out of memory, increase `#SBATCH --mem=8G` in `run_prody.slurm`.
- If a compute node cannot access the internet, keep using `--no-download` and make sure `structures/` exists before submitting.
- To calculate fewer modes for testing, use `--n-modes 20`; for the full VDOS, use `--n-modes all`.
