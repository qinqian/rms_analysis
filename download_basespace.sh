#!/bin/bash -ex

## some case examples: https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-examples

## a case for downloading basecalls
bs download run --id 191551362 --output .

## a case for downloading fastqs, downloading individual fastq.gz would failed in very small files
bs download project  -i 144988846 -o . --extension=fastq.gz

