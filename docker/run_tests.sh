#!/bin/bash

cd /SISRS
pip install -e .
python --version
echo $PWD
ls
pytest -s tests/
