<!--
SPDX-FileCopyrightText:  PyPSA-ASEAN, PyPSA-Earth and PyPSA-Eur Authors

SPDX-License-Identifier: CC-BY-4.0
-->
# Run Tutorials

The tutorial models are tailored for new users to become familiar with PyPSA-ASEAN. They rely on a custom databundle and cutouts that span one week rather than a full year. Along with a reduced number of clusters, this setup allows the model to be run using open-source solvers such as HiGHS.

## Steps to Create and Run Tutorials

To create these tutorial models, follow these steps:

### 1. Activate Your Environment

Activate the conda environment as described in [Installation](installation.md)

```bash
$ conda activate pypsa-earth
```

### 2. Define Your Tutorial Configuration

Predefined configurations are available in `configs/tutorials/`, each representing a different model setup:

* `config.asean.yaml`: A simplified regional model based on `configs/config.asean.yaml`.
* `config.java-bali.yaml`: An Indonesian model focused on the Java–Bali power network.
* `config.peninsular.yaml`: A Malaysian model centered on the Malay Peninsula network, also incorporating selected neighbouring countries.

These are example snippets from one of the tutorial configurations:

```yaml
# configs/tutorials/config.asean.yaml

tutorial: true

run:
  name: baseline-aims-3H-tutorial # use this to keep track of runs with different settings
  sector_name: baseline-aims-3H-tutorial

scenario:
  clusters: [50]
  planning_horizons:
  - 2030
  - 2040
  - 2050

snapshots:
  start: "2013-03-01"
  end: "2013-03-07"
...
```

### 3. Run the Tutorial Configuration

To run this configuration:

```bash
snakemake solve_sector_networks --configfile configs/config.asean.yaml configs/tutorials/config.asean.yaml -call
```

* `solve_sector_networks` specifies the rule to execute
* `--configfile` indicates the configuration file that overrides the default PyPSA-Earth settings. `configs/config.asean.yaml` are used to fill ASEAN specific configuration not defined in the tutorials
* `-call` defines the number of CPUs to use (e.g. `-c8`)

To remove all generated outputs, use the `--delete-all-output` flag:

```bash
snakemake solve_sector_networks --configfile configs/config.asean.yaml configs/tutorials/config.asean.yaml --delete-all-output -call
```

These are extra flags that can help with debugging:

* `-n` performs a dry run
* `--rerun-trigger mtime` reuses existing files and avoids unnecessary reruns
* `--rerun-incomplete` reruns jobs with incomplete or corrupted outputs

### 4. Expand the Model

Once you have successfully run a tutorial scenario, you can expand the model by taking the following steps:

- Configure the model to include the regions you want to analyze.
- Set `tutorial: false` to download the full databundle.
- Increase the spatial and temporal resolution. To extend the temporal resolution, update the configuration as follows:

```yaml
atlite:
  default: asean-2013-era5 # (~7.3 GB) replaces asean-2013-era5-tutorial
```

- (Experimental) Create your own custom adjustments in `final_asean_adjustment.py`
- (Experimental) Enable additional sector coupling by setting:

```yaml
final_adjustment:
  only_elec_network: false
```

- If the solving time is too long, consider using the `gurobi` solver to improve performance by updating the configuration as follows:

```yaml
solving:
  solver:
    name: gurobi
    options: gurobi-default
```

If you encounter issues while building or expanding the model, please report them by opening an issue in the PyPSA-ASEAN GitHub repository.
