---
label: Install arachne
icon: desktop-download
order: 99
---

>>> Clone the arachne repository
```bash
git clone --recursive https://github.com/pdimens/arachne.git
```
The inclusion of `--recursive` is important to make sure `bwa` and `jemalloc` dependencies are cloned as well. 

>>> Execute the makefile
=== Direct compilation
Direct compilation requires a few dependencies in your software environment:
- autoconf
- automake
- c-compiler
- go >= 1.9.2
- zlib

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
