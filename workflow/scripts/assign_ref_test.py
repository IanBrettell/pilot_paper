# Send stdout and stderr to log file
import sys,os
import logging, traceback
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    )
def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(''.join(["Uncaught exception: ",
                         *traceback.format_exception(exc_type, exc_value, exc_traceback)
                         ])
                 )
# Install exception handler
sys.excepthook = handle_exception

# Import libraries

import pandas as pd

# Get variables

## Debug
#IN_FILE = "/hps/nobackup/birney/users/ian/pilot/split/open_field/session_20190611_1331_icab_icab_R_q1/trajectories_wo_gaps/trajectories_wo_gaps.trajectories.csv"
#REF_LOC = "NA"
#FPS = 30
#OUT_FILE = "/nfs/research/birney/users/ian/MIKK_F2_tracking/novel_object/20211118_1110_R_q3.csv"

## True
IN_FILE = snakemake.input[0]
REF_LOC = snakemake.params.ref_loc
FPS = int(snakemake.params.fps)
OUT_FILE = snakemake.output[0]

# Read in .csv

df = pd.read_csv(IN_FILE)

# Remove hash from first column name

df.rename(columns = {'# x1':'x1'}, inplace = True)

# Remove NaNs (so that the first row isn't empty)

df_complete = df.dropna().reset_index()

# Get column name matching REF location

if REF_LOC == "Left" or "Right":
    # Pull out first row
    target = df_complete.loc[0,['x1', 'x2']]
    # Get column name with min or max value based on REF location
    if REF_LOC == "Left":
        ref_col = target.idxmin()
    elif REF_LOC == "Right":
        ref_col = target.idxmax()

if REF_LOC == "Top" or "Bottom":
    # Pull out first row
    target = df_complete.loc[0,['y1', 'y2']]
    # Get column name with min or max value based on REF location
    if REF_LOC == "Top":
        ref_col = target.idxmin()
    elif REF_LOC == "Bottom":
        ref_col = target.idxmax()

# If both fishes are from the same line, make the first object the reference 
if REF_LOC == "NA":
    ref_col = 'x1'

# Rename columns

if ref_col[1] == '1':
    df.rename(columns = {'x1':'ref_x', 'y1':'ref_y', 'x2':'test_x', 'y2':'test_y'}, inplace = True)
elif ref_col[1] == '2':
    df.rename(columns = {'x1':'test_x', 'y1':'test_y', 'x2':'ref_x', 'y2':'ref_y'}, inplace = True)

# Add `frame` and `seconds` column

df_final = df.copy()
df_final['frame'] = [i for i in range(1,len(df) + 1)]
df_final['seconds'] = df_final['frame'] / FPS

# Filter for frames up to 10 minutes

df_final = df_final.loc[(df_final['seconds'] <= 600)]

# Reorder columns

df_final = df_final[['frame', 'seconds', 'ref_x', 'ref_y', 'test_x', 'test_y']]

# Write to file

df_final.to_csv(OUT_FILE, index = False)    
