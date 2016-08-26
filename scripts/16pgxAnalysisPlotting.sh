#!/bin/bash
ROOT_DIR=..
DATA=$ROOT_DIR/data
PGX=$ROOT_DIR/pgx_analyses

mkdir -p $PGX/plots

# Run main analysis in R
Rscript pgx1000gpMain.R 
