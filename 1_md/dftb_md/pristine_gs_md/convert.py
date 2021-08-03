
def dftb_traj2xyz_traj(in_filename, out_filename):    
    """
    This is an auxiliary function to convert the DFTB+-generated extended xyz trajectory file
    to a more coomon xyz trajectory file format

    Args:
        in_filename ( string ):  the name of the input file
        out_filename ( string ):  the name of the putput file

    Example:
        dftb_traj2xyz_traj("md.xyz", "md_reduced.xyz")

    """

    f = open(in_filename, "r")
    A = f.readlines()
    f.close()

    f = open(out_filename, "w")
    for line in A:
        tmp = line.split()
        res = line

        if len(tmp) == 8:
            res = F"  {tmp[0]}  {tmp[1]}  {tmp[2]}  {tmp[3]}\n"

        f.write(res)
    f.close()


dftb_traj2xyz_traj("md.xyz", "md_reduced.xyz")

