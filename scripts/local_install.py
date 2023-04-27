import os
from setuptools import setup

#Kma
cmd = "git clone https://bitbucket.org/genomicepidemiology/kma.git"
os.system(cmd)
cmd = "cd kma && make"
os.system(cmd)
cmd = "cp kma ~/bin/. && cd .."
os.system(cmd)
#ccphylo
cmd = "git clone https://bitbucket.org/genomicepidemiology/ccphylo.git"
os.system(cmd)
cmd = "cd ccphylo && make"
os.system(cmd)
cmd = "cp ccphylo ~/bin/. && cd .."
os.system(cmd)