#!/bin/bash
version=$(cat version)
packagename="qpa-$version"
tarfile="$packagename.tar.gz"
tar --transform "s,^,$packagename/," -zcf $tarfile \
    CHANGES LICENSE README PackageInfo.g init.g read.g version \
    lib/*.{gd,gi,g} examples/*.g tst/*.{g,tst} \
    doc/*.{xml,txt,html,css,js} doc/manual.pdf doc/manual.six doc/MakeQPADoc.gi
