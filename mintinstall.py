import  os
cmd = "chmod a+x Mint3"
os.system(cmd)
cmd = "mv Mint3 ~/bin/"
os.system(cmd)
cmd = "chmod a+x Mint3Functions.py"
os.system(cmd)
cmd = "mv Mint3Functions.py ~/bin~"
os.system(cmd)
#Kma
cmd = "git clone https://bitbucket.org/genomicepidemiology/kma.git"
os.system(cmd)
cmd = "cd kma && make"
os.system(cmd)
cmd = "mv kma ~/bin/"
os.system(cmd)
cmd = "cd .."
os.system(cmd)
#ccphylo
cmd = "git clone https://bitbucket.org/genomicepidemiology/ccphylo.git"
os.system(cmd)
cmd = "cd ccphylo && make"
os.system(cmd)
cmd = "mv ccphylo ~/bin/"
os.system(cmd)
cmd = "cd .."
os.system(cmdos.system(cmd)
