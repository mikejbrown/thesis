#!/bin/sh

if [ -z "$(which lyx | grep 'No lyx in')" ];
then
  lyxPath="C:\\Program Files (x86)\\LyX 2.2\\bin\\lyx.exe"
else
  lyxPath="lyx"
fi

fname=thesis-"$(date '+%d%m%Y-%H%M%S')".pdf

# make sure the version number is up to date
git describe --dirty --always > version

"$lyxPath" -E pdf2 $fname thesis_main.lyx
