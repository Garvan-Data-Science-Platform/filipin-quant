# Running the Filipin GUI on Windows (with a desktop icon)

This folder contains everything needed to launch the GUI by double-clicking an
icon on the Windows desktop.

| File | What it is |
|------|------------|
| `launch_filipin_gui.bat` | The launcher — activates the conda env and starts the GUI. |
| `Create_Desktop_Shortcut.bat` | Double-click **once** to create the icon on **your own** Desktop. |
| `Create_Shortcut_All_Users.bat` | Admin: put the icon on the **shared Public Desktop** (visible to everyone, e.g. deploying for someone else). |
| `Create_Desktop_Shortcut.ps1` | The PowerShell the `.bat` wrappers run (do not run directly). |
| `filipin.ico` | The app icon used by the shortcut. |

## One-time setup

**1. Copy the whole project folder** to the Windows machine (the folder that
contains `filipin_gui.py`, `segment.py`, `filipin_quant.py`, this `windows\`
folder, etc.). Use a short, local, space-free path, e.g. `C:\filipin`. Avoid
OneDrive-synced folders (Desktop/Documents on many managed machines) and network
drives.

**2. Create the conda environment.** The pinned spec is in `environment.yml` in
the project root — from an *Anaconda Prompt* in the project folder:

```bat
conda env create -f environment.yml
conda activate cell_quant_sam
```

This creates an env named **`cell_quant_sam`** with napari + Cellpose(-SAM) + a
GPU torch build. **Read the CUDA note at the top of `environment.yml`** — the
pinned torch is the CUDA 12.8 build; if this machine has an older GPU/driver,
edit that one line to the matching CUDA version before creating the env.

Quick check that the env works (should print `torch cuda: True`):
```bat
conda activate cell_quant_sam
python -c "import napari, cellpose, torch; print('torch cuda:', torch.cuda.is_available())"
```

> If you name the env something other than `cell_quant_sam`, open
> `launch_filipin_gui.bat` in Notepad and change the `ENV_NAME=` line.

**3. Make the desktop icon.** In this `windows\` folder, double-click
**`Create_Desktop_Shortcut.bat`**. A **"Filipin Quant"** icon appears on the
Desktop.

## Installing for another user (admin deploying for a standard user)

Do this when you set the machine up from an admin account but someone else runs
it. **The critical point: where the env is created.**

> By default `conda env create` puts the env under **your admin profile**
> (`C:\Users\admin-...\.conda\envs\`). A standard (non-admin) user **cannot read
> another user's profile**, so an env created there will NOT work for them. You
> must force it into an all-users location.

If the machine already has an all-users conda base (e.g. `C:\ProgramData\miniconda3`),
reinstalling conda is not needed; otherwise install Miniconda "for all users" first. Steps:

1. Put the code in a shared local path, e.g. `C:\filipin` (step 1 above).
2. Open **Anaconda Prompt** → **Run as administrator** (so `C:\ProgramData` is
   writable). Make sure the active base is the system one:
   ```bat
   conda activate C:\ProgramData\miniconda3
   conda info --base        :: should print C:\ProgramData\miniconda3
   ```
3. Create the env **with an explicit all-users prefix** (this is what makes it
   reachable by the end user):
   ```bat
   cd C:\filipin
   conda env create -f environment.yml -p C:\ProgramData\miniconda3\envs\cell_quant_sam
   ```
   (The `-p` path overrides the `name:` in the yml. Do NOT create it by name
   only — that lands in your private `.conda\envs`.)
4. Tell the launcher to use that env by path: open `launch_filipin_gui.bat` in
   Notepad and set the prefix line:
   ```bat
   set "ENV_PREFIX=C:\ProgramData\miniconda3\envs\cell_quant_sam"
   ```
5. For the icon, run **`Create_Shortcut_All_Users.bat`** → right-click →
   **Run as administrator**. It places "Filipin Quant" on the **Public Desktop**,
   so it shows on the end user's Desktop (and everyone's).
6. **Verify as the end user if you can**: log in as them (or have them try), double-click
   the icon. The first segmentation downloads the Cellpose-SAM model weights into
   *their* profile, so **their first run needs internet once** (a one-off
   ~few-hundred-MB download; it can look like a hang — it isn't).

## Daily use

Double-click **Filipin Quant** on the Desktop. A small console window opens
(showing progress), then the napari GUI appears. Follow the
"Using the GUI" section of the main `README.md`.

Leave the console window open while you work — closing it closes the GUI. It also
prints any errors, which is handy if something goes wrong.

## Troubleshooting

- **"Could not find conda automatically."** The launcher probes the usual install
  paths. Open *Anaconda Prompt*, run `where conda`, and either add conda to your
  PATH or paste that folder into the `for %%D in (...)` list in
  `launch_filipin_gui.bat`.
- **Window flashes and closes / error code shown.** The console stays open on
  error — read the message. Usually the env name is wrong or a package is missing
  (`conda activate cell_quant_sam` then re-run the import check in step 2).
- **GUI opens but segmentation is very slow / "Killed".** The GPU is busy or torch
  isn't seeing CUDA. Confirm `torch.cuda.is_available()` is `True`, or untick
  **use GPU** in the GUI to fall back to CPU.
- **Shortcut has a blank icon.** Make sure `filipin.ico` is still in this folder,
  then re-run `Create_Desktop_Shortcut.bat`.
