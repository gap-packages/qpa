#!/bin/bash
version=$(cat version)
packagename="qpa-$version"
tarfile="$packagename.tar.gz"
releasedir="releases/$packagename"
mkdir -p $releasedir
cp -r \
   CHANGES LICENSE README PackageInfo.g init.g read.g version \
   lib examples tst doc \
   $releasedir
cd releases
COPYFILE_DISABLE=1 tar -zcf $tarfile \
    $packagename/CHANGES $packagename/LICENSE $packagename/README \
    $packagename/PackageInfo.g $packagename/init.g $packagename/read.g $packagename/version \
    $packagename/lib/*.{gd,gi} $packagename/examples/*.g $packagename/tst/*.{g,tst} \
    $packagename/doc/*.{xml,txt,html,css,js} $packagename/doc/manual.pdf \
    $packagename/doc/manual.six $packagename/doc/MakeQPADoc.gi
