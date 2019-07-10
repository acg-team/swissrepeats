#!/bin/sh

pandoc  -o manuscript.pdf --latex-engine=xelatex manuscript.md
pandoc  -o manuscript.html --template pandoc-bootstrap-template-master/template.html --css pandoc-bootstrap-template-master/template.css --self-contained --toc --toc-depth 2 manuscript.md
