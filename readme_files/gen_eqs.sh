#!/bin/bash

# Generate equations in png format
for e in es eo gs go
do 
	pdfcrop $e.pdf
	convert -density 400 $e-crop.pdf $e.png
done

