#!/bin/sh
tarfile=qpa/qpa-version-$(cat version).tar
cd ../
tar -cvf $tarfile qpa/CHANGES qpa/LICENSE qpa/PackageInfo.g qpa/README qpa/create-release.sh qpa/doc qpa/examples qpa/init.g qpa/lib qpa/read.g qpa/tst qpa/version 
gzip $tarfile
