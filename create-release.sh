#!/bin/bash
version=$(cat version)
packagename="qpa-$version"
tarfile="$packagename.tar.gz"
releasedir="releases/$packagename"
mkdir -p $releasedir
cp -r * $releasedir
cd releases
tar -zcf $tarfile \
    $packagename/CHANGES $packagename/LICENSE $packagename/README \
    $packagename/PackageInfo.g $packagename/init.g $packagename/read.g $packagename/version \
    $packagename/lib/*.{gd,gi} $packagename/examples/*.g $packagename/tst/*.{g,tst} \
    $packagename/doc/*.{xml,txt,html,css,js} $packagename/doc/manual.pdf \
    $packagename/doc/manual.six $packagename/doc/MakeQPADoc.gi
