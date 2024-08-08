# Snakemake-Flagger
Workflow for running [`Flagger`](https://github.com/mobinasri/flagger/tree/main) without mapping.

While several workflows in WDL exist, they have poor interoperability with our existing workflows.

Based on Flagger's excellent [documentation](https://github.com/mobinasri/flagger/tree/main/docs/flagger).

### Getting Started
```bash
git clone 
```

### Configuration

### Usage
```bash
snakemake -np -c 1 --configfile config/config.yaml
```

### Module
```python
SAMPLE_NAMES = ["sample_1"]
```