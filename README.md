# 2023 M2 Feelpp Sonorhc 
Author: Christophe Prud'homme ![GitHub](https://github.com/prudhomm)
Version: v2

![CI](https://github.com/feelpp/2023-m2-feelpp-sonorhc/workflows/CI/badge.svg)

This repository provides a basic starting point for a Feel++ application including:

- [x] Feel++ applications in C++ to use Feel++ and Feel++ toolboxes in `src`
- [x] Documentation using asciidoc and antora
- [x] Python Feel++ notebooks that can be downloaded from the documentation
- [x] Continuous integration including tests for the C++ applications
- [x] Docker image generation for the project

The documentation for 2023-m2-feelpp-sonorhc is available [here](https://feelpp.github.io/2023-m2-feelpp-sonorhc), and you can build on it for your project by enabling the [github pages](https://docs.github.com/en/pages) for your repository.

## Renaming the project

By default, the project is named `2023-m2-feelpp-sonorhc` if you cloned the repository `feelpp/2023-m2-feelpp-sonorhc`. However, if you used the previous repository as a template, then the project is renamed using the name of the repository using the script `rename.sh` at the initialization of the repository. If the name does not suit you, you can change it again using the script `rename.sh` and providing the new name as an argument.

> **Warning**: The script `rename.sh` will rename the project. However, some URLs might not be set properly if you rename the project yourself. You need to check the following files: `docs/site.yml` and `docs/package.json` and fix the URLs after the rename process is done.

## Updating the project version

The version of the project is defined in the files `CMakeLists.txt`, `docs/antora.yml`, and `docs/package.json`. You need to update with the same version in all files.

## Release process

- [x] Update the version in `CMakeLists.txt`
- [x] Update the version in `docs/antora.yml`
- [x] Commit the changes with the message "Release vx.y.z". At this point, the CI will generate the docker image and push it to Docker Hub.
