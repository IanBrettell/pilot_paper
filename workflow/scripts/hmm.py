# Send stdout and stderr to log file
import sys,os
import logging, traceback
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(a
                    sctime)s %(message)s',
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

# Import variables

## Debug 
IN_FILE = "/hps/nobackup/birney/users/ian/pilot/with_metrics/novel_object/20190616_1717_icab_icab_L/q1/0.2.csv"
N_STATES = int("5")

## True
IN_FILE = snakemake.input[0]
N_STATES = int(snakemake.params.n_states)

