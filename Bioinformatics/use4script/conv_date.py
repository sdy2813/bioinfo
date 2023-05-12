from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from datetime import datetime, timedelta
import argparse
import re

def date_to_decimal(year, month=1, day=1):
    start = datetime(year, 1, 1)
    if month == 1 and day == 1:
        date = datetime(year, 7, 2)  # the middle of the year
    elif day == 1:
        date = datetime(year, month, 15)  # the middle of the month
    else:
        date = datetime(year, month, day)
    end = datetime(year + 1, 1, 1)
    return year + ((date - start).days / (end - start).days)

def process_records(input_file, output_file):
    records = []
    for record in SeqIO.parse(input_file, "fasta"):
        # Extract date from the sequence id
        match = re.search(r'@(\d{4})(?:-(\d{2})(?:-(\d{2}))?)?', record.id)
        if match:
            year, month, day = match.groups()
            year = int(year)
            month = int(month) if month else 1
            day = int(day) if day else 1
            decimal_date = date_to_decimal(year, month, day)

            # Replace the date in the sequence id with the decimal date
            new_id = re.sub(r'@\d{4}(?:-\d{2}(?:-\d{2})?)?', f'@{decimal_date:.2f}', record.id)
            new_record = SeqRecord(record.seq, id=new_id, description='')
            records.append(new_record)
    
    SeqIO.write(records, output_file, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help = False, usage = '\npython3 conv_date.py -i [fa] -o [fa]', description='Process fasta file.')
    parser.add_argument('-i','--input', metavar = '[fa]', help = 'input file，fasta format', type=str, required=True)
    parser.add_argument('-o','--output', metavar = '[fa]', help = 'output file，fasta format', type=str, required=True)
    parser.add_argument('-h', '--help', action = 'help', help = 'help info')
    args = parser.parse_args()

    process_records(args.i, args.o)

