import os
import numpy as np
import pandas as pd
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile


def get_sorted_times(dir):
    """
    Get the sorted time directories in a given directory,
    e.g. postProcessing/wallHeatFlux/
    """
    # return [float(string) for string in os.listdir(dir)]
    times = os.listdir(dir)
    return sorted(times, key=lambda x: float(x))

def read_parameters(path="parameters"):
    """
    Read the parameters file, in OpenFOAM format, using PyFoam.
    """
    return ParsedParameterFile(path, noHeader=True).getValueDict()


def load(filename, start=0, end=-1):
    # data = np.loadtxt(filename, comments="#")
    # times = data[start:end, 0]
    # values = [data[start:end, i] for i in range(1, data.shape[1])]
    # return (times, *values)
    """
    Rewrite using pandas. 
    """
    with open(filename, "r") as f:
        is_comment = True
        prev_line = ""
        while is_comment:
            line = f.readline()
            if line[0] != "#":
                is_comment = False
                column_names = prev_line[1:].split()
            prev_line = line
    df = pd.read_csv(filename, comment="#", names=column_names, sep="\s+")
    return df

def load_all_times(times_dir, use_latest_run=True):
    """
    Load data from all time directories in a given directory.
    - times_dir: Directory containing time directories, e.g. "postProcessing/forces"
    - use_latest_run: If True, only load the latest run in each time directory, e.g. "wallHeatFlux_0.dat"
    Returns a merged pandas DataFrame with all data.
    """
    times = get_sorted_times(times_dir)
    data = pd.DataFrame()

    for time in times:
        if use_latest_run:
            files = os.listdir(f"{times_dir}/{time}")
        else:
            files = os.listdir(f"{times_dir}/{time}")
            files = [f for f in files if not f.endswith(f"_{time}")]
        for file in files:
            filename = f"{times_dir}/{time}/{file}"

            with open(filename, "r") as f:
                is_comment = True
                prev_line = ""
                while is_comment:
                    line = f.readline()
                    if line[0] != "#":
                        is_comment = False
                        column_names = prev_line[1:].split()
                    prev_line = line
            df = pd.read_csv(filename, comment="#", names=column_names, sep="\s+")
            # df["time"] = float(time)
            data = pd.concat([data, df], axis=0)
    # Sort by time
    data = data.sort_values("Time")
    return data


def load_wallHeatFlux(filename, start=0, end=-1):
    wall_heat_flows = {}
    times = set()

    # Find unique patch names
    lines = open(filename, "r").readlines()
    for line in lines:
        if line[0] != "#":
            line = line.rstrip().split("\t")
            patch_name = line[1]
            if patch_name not in wall_heat_flows:
                wall_heat_flows[patch_name] = []
            else:
                break

    # Read wall heat flows
    lines = open(filename, "r").readlines()
    for line in lines[
        : -len(wall_heat_flows)
    ]:  # Last time is for some reason written twice, so don't read last entries
        if line[0] != "#":
            line = line.rstrip().split("\t")
            time = float(line[0])
            patch_name = line[1]
            times.add(time)
            wall_heat_flows[patch_name].append(float(line[-1]))

    times = np.array(sorted(list(times)))[start:end]
    for key in wall_heat_flows:
        wall_heat_flows[key] = np.array(wall_heat_flows[key])[start:end]

    return times, wall_heat_flows


def load_sampling_set(sampling_dir, time, set_name, org=False):
    """
    Load data from OpenFOAM sample sets.
    Parameters:
    - sampling_dir: Name of the sampling directory, e.g. "postProcessing/linesampling"
    - time: Time directory, e.g. "0.1"
    - set_name: Name of the set, e.g. "centerline"

    Returns:
    - pandas DataFrame with the data
    """
    files = os.listdir(f"{sampling_dir}/{time}/")
    samples = [f.split("_")[0] for f in files]
    file_type = ''


    for file in files:
        print(f"Checking if {set_name} is in {file}")
        print(f"{set_name}" in file)

        check = f"{set_name}" in file if org else f"{set_name}_" in file

        if check:
            print(f"Loading {file}")
            filename = f"{sampling_dir}/{time}/{file}"
            filetype = file.split(".")[1]
            
            if filetype == "csv":
                df = pd.read_csv(filename)
                df.columns = [str(col) for col in df.columns]
                if not org:
                    df.columns = [col.replace("_0", "x").replace("_1", "y").replace("_2", "z") for col in df.columns]
                else:
                    df.columns = [col.replace("_x", "x").replace("_y", "y").replace("_z", "z") for col in df.columns]
                print(f"Loaded {filename}")


            elif filetype == "xy":
                fields_string = file.replace(f"{set_name}_", "").split(".")[0]
                fields = ["x"] + fields_string.split("_")
                raw_data = np.loadtxt(filename)
                data = {}

                # Split vector fields into components
                offset = 0
                for i in range(len(fields)):
                    j = i + offset
                    field_name = fields[i]
                    if "U" in field_name:
                        data[f"{field_name}x"] = raw_data[:, j]
                        offset += 1
                        j = i + offset
                        data[f"{field_name}y"] = raw_data[:, j]
                        offset += 1
                        j = i + offset
                        data[f"{field_name}z"] = raw_data[:, j]
                    else:
                        data[field_name] = raw_data[:, j]

                df = pd.DataFrame(data)

            # Remove rows that have duplicated values in the first column
            df.drop_duplicates(subset=df.columns[0], inplace=True)

            return df

    return data

