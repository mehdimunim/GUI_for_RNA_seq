#!/bin/bash
inputFile=../QM6aAnnotationIFPEN2021strict.gff
echo "format "$inputFile
sed -i 's/|/_/' $inputFile  