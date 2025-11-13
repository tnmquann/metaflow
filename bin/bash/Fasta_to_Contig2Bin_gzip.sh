#!/usr/bin/env bash

function display_help() {
    echo " "
    echo "Fasta_to_Scaffolds2Bin: Converts genome bins in fasta format (including .gz) to scaffolds-to-bin table."
    echo " (DAS Tool helper script)"
    echo " "
    echo "Usage: Fasta_to_Scaffolds2Bin.sh -e fasta > my_scaffolds2bin.tsv"
    echo " "
    echo "   -e, --extension            Extension of fasta files. (default: fasta)"
    echo "   -i, --input_folder         Folder with bins in fasta format. (default: ./)"
    echo "   -h, --help                 Show this message."
    echo " "
    exit 1
}

extension="fasta"
folder="."

while [ "$1" != "" ]; do
    case $1 in
        -e | --extension )      shift
                                extension=$1
                                ;;
        -i | --input_folder )   shift
                                folder=$1
                                ;;
        -h | --help )           display_help
                                exit
                                ;;
        * )                     display_help
                                exit 1
    esac
    shift
done

# Loop through all fasta and gzipped fasta files
for i in "$folder"/*."$extension"*; do
    [ -e "$i" ] || continue  # Skip if no matching files found

    binname=$(basename "$i" | sed -E "s/\\.$extension(\\.gz)?//g")

    # Detect gzipped vs plain file
    if [[ "$i" == *.gz ]]; then
        gzip -cd "$i" | grep ">" | cut -f1 | perl -pe "s/>//g" | perl -pe "s/\n/\t$binname\n/g"
    else
        grep ">" "$i" | cut -f1 | perl -pe "s/>//g" | perl -pe "s/\n/\t$binname\n/g"
    fi
done
