import  os
cmd = "chmod a+x MinIONTree"
os.system(cmd)
cmd = "mv MinIONTree ~/bin/"
os.system(cmd)
cmd = "chmod a+x MinIONTreeFunctions.py"
os.system(cmd)
cmd = "mv MinIONTreeFunctions.py ~/bin/"
os.system(cmd)
#Kma
cmd = "git clone https://bitbucket.org/genomicepidemiology/kma.git"
os.system(cmd)
cmd = "cd kma && make"
os.system(cmd)
cmd = "mv kma/kma ~/bin/"
os.system(cmd)
cmd = "cd .."
os.system(cmd)
#ccphylo
cmd = "git clone https://bitbucket.org/genomicepidemiology/ccphylo.git"
os.system(cmd)
cmd = "cd ccphylo && make"
os.system(cmd)
cmd = "mv ccphylo/ccphylo ~/bin/"
os.system(cmd)
cmd = "cd .."
os.system(cmd)
