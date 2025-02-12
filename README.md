# GenEvol

GenEvol is a likelihood-based method for estimating gene family dynamics, including gene gain, loss, innovation, and duplication in plants. It provides a robust framework for analyzing evolutionary changes in gene families across different plant species.

## Installation

To install GenEvol, run the following script:

```sh
script/install_sources.sh --workdir <your_workdir>
```

If you're using codespaces, you can omit the workdir flag and use the default `/workspaces` directory.
Once the sources are installed, execute the following commands:

```sh
script/build.sh --workdir <your_workdir>
```