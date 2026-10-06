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
The inclusion of `--recursive` is important to make sure the `minibwa` dependency is cloned as well. 

>>> Execute the makefile
=== Direct compilation
Direct compilation requires a few dependencies in your software environment:
- Go
- a C compiler and zlib
- a CPU with SSE4.2 (x86-64) or NEON (arm64), which minibwa requires
- optionally OpenMP (`libgomp`), used by `arachne index` for multi-threaded index construction. It is detected
  automatically when building; without it indexing is single-threaded

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

This results in the compiled executable binaries `bin/arachne` and `bin/minibwa`. You can use them there
or copy them into another path.
>>>
