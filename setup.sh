#!/usr/bin/env bash

set -ex

PWD="$(dirname "$(realpath "$0")")"

cd "${PWD}"
npm install
python setup.py install --user

#python -m eel main.py static