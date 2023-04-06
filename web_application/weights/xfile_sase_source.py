# XFILE SASE SOURCE
#   Takes SASE temporal profile and creates a CRETIN source 
#   xfile that can be appended to a generator using the #include command.
# INPUT:
#   HISTID unique history id for the source command.
#   FILENAME name of the input SASE file.
#   TIMESHIFT shift the profile time by specified amount.
# OUTPUT:
#   xfile with source commad and specified history id.
# USAGE:
# for file in 9fs*.dat; do python xfile_sase_source.py $file ; done

import sys
from matplotlib import pyplot


def cretin_dump_profile(prefix: str, histid: int, time: list, power: list):
    outfile = open(f"{prefix}.xfile", 'w')
    source = f"source jbndry 1 E1 E2 value history {histid} MULT\n"
    pulse = f"history {histid} ILASER 1e-15\n"
    outfile.writelines([source, pulse])
    for index in range(1, len(time)):
        tvline = f"tv {time[index-1]:.4f} {power[index-1]:.4f}\n"
        outfile.write(tvline)
    outfile.close()
    return


def read_bin_profile(file_name: str, x_translation: float):
    infile = open(file_name, 'r')
    profile = infile.read()
    infile.close()
    time = []
    power = []
    for line in profile.split("\n"):
        if line.startswith("#") or len(line) == 0:
            continue
        split_line = line.split()
        time.append(float(split_line[0])+x_translation)
        power.append(float(split_line[1])*1E-9)

    # normalize power
    power = [p/max(power) for p in power]
    return time, power


def main(file_name: str, histid: int, x_translation: float, show_plot: bool, write_file: bool):
    time, power = read_bin_profile(file_name, x_translation)
    print(f"timestep:{round(time[1]-time[0], 4)} [fs]")
    prefix = file_name.split(".")[0]

    if show_plot:
        pyplot.figure()
        pyplot.subplot(1, 1, 1)
        pyplot.plot(time, power, color="tab:orange")
        pyplot.ylabel("power [GW]")
        pyplot.xlabel("time [fs]")
        pyplot.show()
    if write_file:
        cretin_dump_profile(prefix, histid, time, power)


if __name__ == "__main__":
    FILE_NAME = '/Users/harry/Downloads/readme.dat'
    TIMESHIFT = 15.0  # fs
    SHOW_PLOT = True
    WRITE_FILE = False
    HISTORYID = 1
    main(FILE_NAME, HISTORYID, TIMESHIFT, SHOW_PLOT, WRITE_FILE)
