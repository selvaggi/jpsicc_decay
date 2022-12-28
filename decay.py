import numpy as np
import math
import random
import sys
import os
from src.BoostDecay import Decay

# ___________________________________________________________________________________________________
def mass(p):
    m2 = p[0] * p[0] - p[1] * p[1] - p[2] * p[2] - p[3] * p[3]
    if m2 > 0:
        return math.sqrt(m2)
    else:
        0.0


# ___________________________________________________________________________________________________
def part_line(
    pid,
    status,
    parent1,
    parent2,
    color1,
    color2,
    px,
    py,
    pz,
    e,
    m,
    v,
    s,
):
    line = "{:>9}  {}    {}    {} {:>4} {:>4} {:+.10e} {:+.10e} {:+.10e} {:.10e} {:.10e} {:.4e} {:.4e} \n".format(
        pid, status, parent1, parent2, color1, color2, px, py, pz, e, m, v, s
    )
    return line


# __________________________________________________________________________________________________
def main():

    # input LHE file (undecayed)
    inputFile = sys.argv[1]

    # output LHE file (decayed)
    outputFile = sys.argv[2]

    # option (h or z decay supported)
    option = sys.argv[3]

    # option (h or z decay supported)
    weight_dir = os.path.abspath(sys.argv[4])

    if len(sys.argv) < 5:
        print(" Usage: decay.py [INPUTLHE] [INPUTLHE] [Z or H] [EVTS/WEIGHTS DIR]")
        sys.exit(1)

    if option != "Z" and option != "H":
        print(" only Z or H decay option are supported")
        sys.exit(1)

    pid = 25
    decay = 0
    SDCfactor = 1

    if option == "H":
        SDCfactor = 1.0994
    elif option == "Z":
        pid = 23
        SDCfactor = 0.26462

    # instantiate decay class
    decay = Decay(weight_dir, SDCfactor)

    spins = [1, -1]

    # open the input LHE file for reading
    with open(inputFile, "r") as infile:

        # Read the file and store the lines in a list
        lines = infile.readlines()

        # open the output file for writing
        for i, line in enumerate(lines):

            # check if the higgs
            if line.startswith("       {}  1    1    2".format(pid)):

                ### access Higgs kinemtatics and info
                string_list = line.split()
                pid = string_list[0]
                status = string_list[1]
                parent1 = string_list[2]
                parent2 = string_list[3]
                color1 = string_list[4]
                color2 = string_list[5]
                px = float(string_list[6])
                py = float(string_list[7])
                pz = float(string_list[8])
                e = float(string_list[9])
                m = float(string_list[10])
                v = string_list[11]
                s = string_list[12]

                p4 = np.array([e, px, py, pz])

                ## compute decay info
                decay_info = decay.Decay_Events(p4)
                c_p4 = decay_info[1]
                cbar_p4 = decay_info[2]
                mum_p4 = decay_info[3]
                mup_p4 = decay_info[4]
                weight = decay_info[5]

                random.shuffle(spins)
                spin_c = spins[0]
                spin_mu = spins[1]

                ## access event info, replace number of particles, event weight and infer mother position
                pos = -1
                for j in range(i):
                    if not "<event>" in lines[i - j]:
                        pos += 1
                    else:
                        break

                # Get the next line and print it
                event_info = lines[i - pos]
                event_info_list = event_info.split()
                number_of_particles = int(event_info_list[0])
                number_of_particles += 4

                event_info_list[0] = "{}      ".format(number_of_particles)
                event_info_list[2] = "{:+.7e}".format(weight)

                new_event_info_line = " ".join(event_info_list)
                new_event_info_line = "{}\n".format(new_event_info_line)

                lines[i - pos] = new_event_info_line

                ###
                mother = part_line(
                    pid,
                    2,
                    1,
                    2,
                    0,
                    0,
                    p4[1],
                    p4[2],
                    p4[3],
                    p4[0],
                    mass(p4),
                    0.0,
                    0,
                )

                c = part_line(
                    4,
                    1,
                    pos,
                    pos,
                    503,
                    0,
                    c_p4[1],
                    c_p4[2],
                    c_p4[3],
                    c_p4[0],
                    mass(c_p4),
                    0.0,
                    spin_c,
                )

                cbar = part_line(
                    -4,
                    1,
                    pos,
                    pos,
                    0,
                    503,
                    cbar_p4[1],
                    cbar_p4[2],
                    cbar_p4[3],
                    cbar_p4[0],
                    mass(cbar_p4),
                    0.0,
                    spin_c,
                )

                mum = part_line(
                    13,
                    1,
                    pos,
                    pos,
                    0,
                    0,
                    mum_p4[1],
                    mum_p4[2],
                    mum_p4[3],
                    mum_p4[0],
                    mass(mum_p4),
                    0.0,
                    spin_mu,
                )

                mup = part_line(
                    -13,
                    1,
                    pos,
                    pos,
                    0,
                    0,
                    mup_p4[1],
                    mup_p4[2],
                    mup_p4[3],
                    mup_p4[0],
                    mass(mup_p4),
                    0.0,
                    spin_mu,
                )

                lines[i] = mother
                lines.insert(i + 1, c)
                lines.insert(i + 2, cbar)
                lines.insert(i + 3, mum)
                lines.insert(i + 4, mup)

    with open(outputFile, "w") as outfile:
        outfile.writelines(lines)


# __________________________________________________________________________________________________
if __name__ == "__main__":
    main()
