#!/bin/bash

# NOTE: required for click library to function with Python 3. May need to use
# a different CLI system. See: http://click.pocoo.org/5/python3/
export LC_ALL=C.UTF-8
export LANG=C.UTF-8

cd /SISRS
pip install -e .
pytest -s tests/
