import csv
import itertools
from pyteomics import mass, parser, fasta

peptide_iterator = itertools.chain.from_iterable(
    parser.cleave(protein.sequence, 'trypsin')
    for protein in fasta.FASTA("./yeast_proteins.fa")
)

records = []

def filter_peptide(x: str):
    return "X" not in x and "B" not in x and "Z" not in x and len(x) >= 5

# Add each peptide's product ions to the index, keyed to the peptide's
# ascending mass order.
for j, peptide in enumerate(
        sorted(filter(filter_peptide, peptide_iterator), key=mass.fast_mass)):
    if len(peptide) < 3:
        continue
    peptide_mass = mass.fast_mass(peptide)
    record = {"type": "PEPTIDE", "sequence": peptide, "mass": peptide_mass, "id": j}
    records.append(record)
    for i in range(1, len(peptide)):
        m = mass.fast_mass(peptide[:i], ion_type='b')
        record = {
            "type": "FRAGMENT",
            "sequence": f"b{i}",
            "mass": m,
            "id": j,
        }
        records.append(record)

        m = mass.fast_mass(peptide[i:], ion_type='y')
        record = {
            "type": "FRAGMENT",
            "sequence": f"y{len(peptide) - i}",
            "mass": m,
            "id": j,
        }
        records.append(record)

with open("test_data.csv", 'wt', newline='') as fh:
    writer = csv.DictWriter(fh, fieldnames=['type', 'sequence', 'mass', 'id'])
    writer.writeheader()
    writer.writerows(records)
