import  os
cmd = "chmod a+x MINTGraft"
os.system(cmd)
cmd = "mv MINTGraft ~/bin/"
os.system(cmd)
cmd = "chmod a+x MINTGraftFunctions.py"
os.system(cmd)
cmd = "mv MINTGraftFunctions.py ~/bin/"
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
