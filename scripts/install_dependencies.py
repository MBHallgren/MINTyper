import os
from setuptools import setup
from pathlib import Path

home = str(Path.home())
#Kma
cmd = "git clone https://bitbucket.org/genomicepidemiology/kma.git"
os.system(cmd)
cmd = "cd kma && make"
os.system(cmd)
cmd = "cp kma/kma {}/bin/kma && cd ..".format(home)
os.system(cmd)
#ccphylo
cmd = "git clone https://bitbucket.org/genomicepidemiology/ccphylo.git"
os.system(cmd)
cmd = "cd ccphylo && make"
os.system(cmd)
cmd = "cp ccphylo/ccphylo {}/bin/ccphylo && cd ..".format(home)
os.system(cmd)