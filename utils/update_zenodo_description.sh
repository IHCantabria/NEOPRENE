#!/bin/bash

pandoc -f markdown -t html ../README.md -o README.html
sed 's/\"/\\\"/g' README.html | tr --delete '\n' > README.json
rm README.html
cat Zenodo_File_Template > zenodo.json
printf \" >> zenodo.json
cat README.json >> zenodo.json
printf '\"\n}\n' >> zenodo.json
rm README.json