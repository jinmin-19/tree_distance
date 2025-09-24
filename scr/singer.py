import argparse
import subprocess
import numpy as np
import tskit
import os


dir="/Users/jinminli/research/tree_distance/simulations/"
names=["Ne2e4_100_1e6","Ne2e4_1000_1e6","Ne2e4_3000_1e6","ooa_300_1e6","ooa_3000_1e6"]

for name in names:
    file_prefix=dir+name+"/"+name
    #convert to vcf
    subprocess.run([
        "bcftools", "view",
        "-Ov",
        "-o", file_prefix+".vcf",
        file_prefix+".vcf.gz"
    ],check=True)

    #run_singer for vcf
    subprocess.run(["python",
        "/Users/jinminli/research/tree_distance/SINGER/run_singer.py",
        "-Ne","2e4",
        "-m","1.25e-8",
        "-start","0",
        "-end","1e6",
        "-length","1e6",
        "-vcf_header",file_prefix,
        "-arg_save",f"{dir}{name}/singer_original/singer_{name}_polar099_",
        "-ts_save",f"{dir}{name}/singer_ts/{name}_",
    ],check=True)