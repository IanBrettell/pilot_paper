# Copy data to local

## Activate env

using Pkg
Pkg.activate("envs/julia")
Pkg.instantiate()

## Create data directory on local
if !isdir(joinpath(homedir(), "Documents/Data/20210110_pilot_test_data")
    mkdir("/Users/brettell/Documents/Data/20210110_pilot_test_data")
end

## Copy over documents
```bash
rclone copy \
    google_drive_ebi:pilot_paper/idtrackerai/final \
    ~/Documents/Data/20210110_pilot_test_data
```

## Import into list

# Load DataFrames
using DataFrames, CSV

data_frames = []
for file in readdir("/Users/brettell/Documents/Data/20210110_pilot_test_data", join = true)
    println(file)
end


target_file = readdir("/Users/brettell/Documents/Data/20210110_pilot_test_data", join = true)[1]

# Read in test DF
test_df = DataFrame(CSV.File(target_file))

# Get column names
plot(test_df[1],)
