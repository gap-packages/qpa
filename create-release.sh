#!/bin/sh
tarfile=qpa-version-$(cat version).tar
git archive --prefix=qpa/ HEAD --format=tar > $tarfile
gzip $tarfile 
