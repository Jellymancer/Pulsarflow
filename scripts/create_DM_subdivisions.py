import os
import numpy as np
import argparse
from typing import List
from math import floor


def create_directory_with_file(root_path: str, dm_values: List[float]) -> None:

    # Convert values to string with new lines
    content = '\n'.join(map(str, dm_values))

    # Determine the directory name based on min and max values
    min_val = floor(min(dm_values))
    max_val = floor(max(dm_values))
    file_name = f"DM_{min_val}_{max_val}.txt"

    # Ensure the root directory exists
    os.makedirs(root_path, exist_ok=True)

    # Path for the new file within the new directory
    file_path = os.path.join(root_path, file_name)

    # Write the values to the file
    with open(file_path, 'w') as file:
        file.write(content)


parser = argparse.ArgumentParser(description='Create DM subdivisions for' 
                                 'Peasoup search')
parser.add_argument('-r', '--root_path', help='Root path where the new '
                    'directories are to be created', type=str, default=os.getcwd())
parser.add_argument('-ds', '--dm_start', help='DM start value', type=float,
                    default=0.0)
parser.add_argument('-de', '--dm_end', help='DM max value', type=float,
                    required=True)
parser.add_argument('-n', '--num_splits', help='Number of splits', type=int,
                    default=16)
parser.add_argument('-s', '--step', help='Step size for DM values', type=float,
                    default=0.1)
parser.add_argument('-d', '--decimals', type=int, default=1,
                    help='No. of decimnals to print for the DM values')


if __name__ == "__main__":
    args = parser.parse_args()

    dm_values = np.arange(args.dm_start, args.dm_end+args.step, args.step)
    full_values = np.array([round(x, args.decimals) for x in dm_values])
    chunksize = len(full_values)/args.num_splits
    vallist = np.array_split(full_values, len(full_values) // chunksize)

    for vals in vallist:
        create_directory_with_file(args.root_path, vals)
