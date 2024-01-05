#!/bin/bash

# Function to display help message
show_help() {
    echo "Usage: $0 --query <query>"
    echo "Options:"
    echo "  --help            Display this help message."
    echo "  --query <query>    Specify the query string for searching."
}

# Parse command line options
while [[ $# -gt 0 ]]; do
    case "$1" in
        --help)
            show_help
            exit 0
            ;;
        --query)
            QUERY="$2"
            shift # past argument
            shift # past value
            ;;
        *)    # unknown option
            show_help
            exit 1
            ;;
    esac
done

# Check if query is provided
if [ -z "$QUERY" ]; then
    echo "Error: No query provided."
    show_help
    exit 1
fi

# Check if esearch, elink, efetch are installed
for cmd in esearch elink efetch; do
    if ! command -v $cmd &> /dev/null; then
        echo "Command $cmd not found. Please install entrez-direct."
        exit 1
    fi
done

# Check if kingfisher is installed
if ! command -v kingfisher &> /dev/null; then
    echo "kingfisher not found. Please install kingfisher."
    exit 1
fi

# Check if SRP.acc file exists
if [ ! -f SRP.acc ]; then
    # If not, run commands to create SRP.acc
    esearch -db gds -query "$QUERY" | elink -target sra | efetch -format docsum | grep "SRP" | awk -F "\"" '{print $2}' | sort | uniq > SRP.acc
fi

# Check if Run.list file exists
if [ ! -f Run.list ]; then
    # If not, run commands to create Run.list
    esearch -db sra -query $(cat SRP.acc) | efetch -format runinfo | cut -f1 -d ',' | grep -v "^Run" > Run.list
fi

# Iterate over each entry in Run.list and execute kingfisher get command
while read -r line; do
    time kingfisher get -r "${line}" -m aws-http
done < Run.list

