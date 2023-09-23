# Linear systems with UMFPACK and matrix market

This code solves linear systems with UMFPACK solver and matrices from Matrix Market.

You may install UMFPACK locally or use a Docker container as explained in the next section.

**Note:** Tested on Ubuntu 23.04.

## Local Installation

**1 Install UMFPACK**

```bash
sudo apt-get install libsuitesparse-dev
```

**2 Download some Matrix Market Files**

```bash
bash scripts/download-from-matrix-market.bash
```

**3 Run Example**

```bash
bash ./all.bash
```

## Using Docker

You may use an existent image from `cpmech` or build the Docker image yourself (see script `docker-build-image.bash`).

When using VS Code, the extension [ms-vscode-remote.remote-containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) is recommended.

**1 Reopen this folder in a container**

![VS Code Remote Development](remote-dev-with-vscode.gif)

**2 Download some Matrix Market Files**

```bash
bash scripts/download-from-matrix-market.bash
```

**3 Run Example**

```bash
bash ./all.bash
```
