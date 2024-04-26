import numpy as np
import pandas as pd
from sqlalchemy import create_engine
import sys
import argparse
import os
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Postgres username, password, and database name
POSTGRES_ADDRESS = os.getenv("DB_HOST")  # Insert your DB address if it's not on Panoply
POSTGRES_PORT = os.getenv("DB_PORT")
POSTGRES_USERNAME = os.getenv("DB_USERNAME")
POSTGRES_PASSWORD = os.getenv("DB_PASSWORD")  # Change this to your Panoply/Postgres password
POSTGRES_DBNAME = 'trapum_web'  # Database name
# A long string that contains the necessary Postgres login information
postgres_str = ('mysql+pymysql://{username}:{password}@{ipaddress}:{port}/{dbname}'
                .format(username=POSTGRES_USERNAME,
                        password=POSTGRES_PASSWORD,
                        ipaddress=POSTGRES_ADDRESS,
                        port=POSTGRES_PORT,
                        dbname=POSTGRES_DBNAME))
# Create the connection
cnx = create_engine(postgres_str)

parser = argparse.ArgumentParser(description='Extract Raw fil Files for a given Source name and beam number')
parser.add_argument('-b', '--beam', help='beam name. for eg cfbf00000', type=str)
parser.add_argument('-t', '--source_name', help='MGPS Source Name', default='MSGPS_L_0243', type=str)
parser.add_argument('-us', '--utc_start', help='UTC start date (YYYY-MM-DD)', type=str, default=None)
parser.add_argument('-ue', '--utc_end', help='UTC end date (YYYY-MM-DD)', type=str, default=None)

args = parser.parse_args()

source_name = args.source_name
file_type_pattern = 'filterbank-raw%'
utc_start = args.utc_start
utc_end = args.utc_end

def sort_and_join(files):
    return ' '.join(sorted(files.split()))

date_condition = ""
if utc_start and utc_end:
    date_condition = "AND p.utc_start BETWEEN %s AND %s "

if args.beam:
    beam_name = args.beam
    query_get_filterbank_files = '''
    SELECT DISTINCT p.id AS pointing_id, t.source_name as target, b.id AS beam_id, b.name AS beam_num, p.utc_start, 
    concat(dp.filepath, "/", dp.filename) AS filterbank_files, ft.name AS file_type
    FROM data_product dp
    JOIN pointing p ON dp.pointing_id = p.id
    JOIN beam b ON dp.beam_id = b.id
    JOIN target t ON p.target_id = t.id
    JOIN file_type ft ON dp.file_type_id = ft.id
    WHERE dp.available = 1 AND t.source_name LIKE %s AND b.name LIKE %s AND ft.name LIKE %s ''' + date_condition + '''order by p.utc_start desc
    '''

    params = (source_name, beam_name, file_type_pattern)
    if utc_start and utc_end:
        params += (utc_start, utc_end)

else:
    query_get_filterbank_files = '''
    SELECT DISTINCT p.id AS pointing_id, t.source_name as target, b.id AS beam_id, b.name AS beam_num, p.utc_start, 
    concat(dp.filepath, "/", dp.filename) AS filterbank_files, ft.name AS file_type
    FROM data_product dp
    JOIN pointing p ON dp.pointing_id = p.id
    JOIN beam b ON dp.beam_id = b.id
    JOIN target t ON p.target_id = t.id
    JOIN file_type ft ON dp.file_type_id = ft.id
    WHERE dp.available = 1 AND t.source_name LIKE %s AND ft.name LIKE %s ''' + date_condition + '''order by p.utc_start desc, b.name asc
    '''

    params = (source_name, file_type_pattern)
    if utc_start and utc_end:
        params += (utc_start, utc_end)

printable_query = query_get_filterbank_files % tuple(repr(p) for p in params)
print(printable_query)

df = pd.read_sql_query(query_get_filterbank_files, con=cnx, params=params)

df.to_csv('filterbank_files.csv', index=False)
df['filterbank_files'] = df['filterbank_files'].apply(lambda x: ' '.join(sorted(x.split())))

# Group by 'pointing_id' and 'beam_id', then sort and concatenate file paths
# Group by 'pointing_id' and 'beam_id', and aggregate
grouped = df.groupby(['pointing_id', 'beam_id']).agg({
    'filterbank_files': lambda x: sort_and_join(' '.join(x)),
    'target': 'first',
    'beam_num': 'first',
    'utc_start': 'first'
}).reset_index()

grouped.to_csv('grouped_files.txt', index=False, header=True, sep='\t')

