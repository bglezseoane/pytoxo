#!/usr/bin/env sh

## PyToxo use example as command line interface
## 
## In this simple shell script we illustrate the use of PyToxo, as command line 
## interface (CLI), with some of the models saved within this repository.
## 
## Warning: run this script from the directory that contains it.


## EDIT HERE
## Comment next line to use the installed PyToxo, uncomment it to use the 
## repository's sources. Useful to be able to use this script during 
## development.
alias pytoxo="python -m pytoxo_cli"


cd ..  # Set path as repository home to access models easier

# Help menu
echo "Printing the program help menu:"
pytoxo --help

# Example 1
echo
echo "Example 1: 'pytoxo models/additive_3.csv 0.2 0.4 0.4 0.4 --max_her':"
pytoxo models/additive_3.csv 0.2 0.4 0.4 0.4 --max_her

# Example 2
tmp_file=$(mktemp)  # Create temporary file to save a table
echo
echo 'Example 2, saving the table to a temporary file: '\''pytoxo models/additive_2.csv --max_her 0.55 0.1 0.2 > "$tmp_file"'\'':'
pytoxo models/additive_2.csv --max_her 0.55 0.1 0.2 > "$tmp_file"
echo 'Printing the content of the temporary file with: '\''cat "$tmp_file"'\'':'
cat "$tmp_file"
rm "$tmp_file"  # Clean

# Example 3
echo
echo "Example 3: 'pytoxo models/threshold_5.csv 0.98 --max_prev 0.1 0.2 0.3 0.25 0.4':"
pytoxo models/threshold_5.csv 0.98 --max_prev 0.1 0.2 0.3 0.25 0.4

# Example 4
echo
echo "Example 4, using GAMETES format: 'pytoxo models/threshold_5.csv 0.98 --max_prev 0.1 0.2 0.3 0.25 0.4 --gametes':"
pytoxo models/threshold_5.csv 0.98 --max_prev 0.1 0.2 0.3 0.25 0.4 --gametes
