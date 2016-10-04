#!/bin/sh
# Prepares an arxiv distributable latex tarball in $buildTarget

if [ -z "$(which lyx | grep 'No lyx in')" ];
then
  lyxPath="C:\\Program Files (x86)\\LyX 2.2\\bin\\lyx.exe"
else
  lyxPath="lyx"
fi

buildTarget="thesis-bundle.tar.gz"
echo
echo "Preparing arXiv ready distribution as: $buildTarget"

# Clean target dir
echo
echo "Cleaning dist/ directory..."
rm -rf dist/
mkdir dist

echo "Exporting lyx -> pdflatex..."
"$lyxPath" -E pdflatex dist/thesis_main.tex thesis_main.lyx

echo
echo "Copying necessary files to dist/..."
cp jcuthesis.cls dist
cp jculogo.pdf dist
mkdir -p dist/code/with-output
cp code/with-output/*.pdf dist/code/with-output

cd dist

# arXiv.org wants \pdfoutput=1 in the latex preamble, but LyX doesn't put it
# there by default. So put it in manually by applying a patch.
cat > fixup.patch <<HERE
3a4
> \pdfoutput=1
HERE
echo
echo "Applying \\pdfoutput=1 patch for arXiv build system..."
patch thesis_main.tex fixup.patch
rm fixup.patch

echo
echo "Building pdf and aux files..."
pdflatex -quiet thesis_main.tex

echo
echo "Building bbl file..."
bibtex -quiet thesis_main

echo
echo "Cleaning auxiliary files..."
rm *.aux *.log *.{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15} *.mp *.t{1,2,3,4,5,6}
rm thesis_main.blg thesis_main.lof thesis_main.lot thesis_main.out thesis_main.toc

echo
echo "Tarballing..."
tar -czvf "$buildTarget" *
mv "$buildTarget" ../"$buildTarget"
cd ..

echo
echo "Cleaning dist/ directory..."
rm -rf dist/

echo
echo "*** All done! ***"
