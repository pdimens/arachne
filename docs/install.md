---
label: Install arachne
icon: desktop-download
order: 99
---

## With `conda`:
This assumes an environment has already been created with `conda create`
> for `mamba`: just swap `conda` with `mamba`
```bash
conda install -c bioconda -c conda-forge bioconda::arachne
```

## With `pixi`
This assumes a pixi project was already created with `pixi init` and has the `conda-forge` and `bioconda` channels added
```bash
pixi add arachne
```

## Manually Compile
>>> Clone the arachne repository
```bash
git clone --recursive https://github.com/pdimens/arachne.git
```
The inclusion of `--recursive` is important to make sure the `bwa` dependency is cloned as well. 

>>> Execute the makefile
=== Direct compilation
Direct compilation requires a few dependencies in your software environment:
- Go
- jemalloc

```bash
cd arachne
make clean; make
```
==- Build with pixi (alternative)
For development portability, Arachne also provides a pixi environment with the build dependencies.
Once arachne is compiled, the pixi environment is no longer needed. This approach assumes pixi is
installed in your software environment. To use the pixi approach:
```bash
cd arachne
pixi run build
```
===

>>>

The compiled and executable `arachne` and `bwa` binaries are now in `bin/`, you can use them there
or copy them into another path. A conda installation will be made available
when the project matures to a stable release.
