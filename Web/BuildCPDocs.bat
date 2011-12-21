@echo off
make latex
cd _build/latex
pdflatex CoolPropdoc.tex
pdflatex CoolPropdoc.tex
pdflatex CoolPropdoc.tex
pdflatex CoolPropdoc.tex
pdflatex CoolPropdoc.tex
copy /Y CoolPropdoc.pdf ..\..\_static\
cd ..\..
make html
