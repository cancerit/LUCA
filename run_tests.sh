#!/bin/sh

PKG=src/luca

pycodestyle --ignore=E501,W504,E266 "${PKG}" tests

pytest --cov="${PKG}" tests/
