---
label: Install arachne
icon: desktop-download
order: 99
---

> [!NOTE]
>Make sure you have a working Go installation (version >= 1.9.2). `go version` should return something like `go version go1.9.2 linux/amd64`.

Clone the arachne repository
```bash
git clone --recursive https://github.com/pdimens/arachne.git
```

Go into the locally cloned repository and execute the makefile
```bash
cd arachne
make
```

The compiled and executable `arachne` and `bwa` binaries are now in `bin/`, you can use them there or copy them into another path.
A conda installation will be made available when the project mature's to a stable release.
