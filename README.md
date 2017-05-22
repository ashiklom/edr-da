# ED and EDR test suite

Author: Alexey Shiklomanov

The `redr` package contains functions to facilitate setting up and running EDR test cases and analyses.
The R scripts in other directories contain examples of different EDR-related analyses.

## Setup

Copy the `config.example.R` to `config.R` and modify the paths as necessary.
The only paths you should _have_ to change for things to work are the executable paths, as the other paths are relative to this repository, and the default should work out of the box.
