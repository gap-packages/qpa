#!/bin/sh
tarfile=qpa-version-$(cat version).tar
git archive HEAD --format=tar > $tarfile
tar -r -f $tarfile doc/chap0.html doc/manual.pdf doc/manual.six 
gzip $tarfile 