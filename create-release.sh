#!/bin/sh
tarfile=qpa-version-$(cat version).tar
git archive HEAD --format=tar > $tarfile
gzip $tarfile 
