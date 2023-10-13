#get full path of installed directory and save somewhere internally
import os

here = os.getcwd()

with open(".proj_dir_path.dat", "w+") as proj_dir:
    proj_dir.write(here)

