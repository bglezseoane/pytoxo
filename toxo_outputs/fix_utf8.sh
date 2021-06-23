#!/usr/bin/env bash

## PyToxo
##
## A Python library for calculating penetrance tables of any
## bivariate epistasis model.
##
## This script convert all Toxo output files of this folder to UTF-8, since
## Toxo does not always encode its outputs respecting this format and for
## PyToxo it is necessary.
##
## Copyright 2021 Borja González Seoane
##

find . -name '*.csv' -exec sh -c "iconv -t UTF-8 {} > {}.utf8"  \; -exec mv "{}".utf8 "{}" \;
