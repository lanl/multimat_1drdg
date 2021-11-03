import os
import subprocess

print("In " + os.getcwd())

# clean directory
fileRange = ["absl1errors.dat", "absl2errors.dat", "logl1errors.dat", "logl2errors.dat"]
for fi in fileRange:
  if os.path.exists(fi):
    os.remove(fi)

# read control file
cntlSrc = open("setflow.cntl","r")
cntlTxt = cntlSrc.read()
cntlSrc.close()

# modify file
cntlLines = cntlTxt.split("\n")
elemRange = [25, 50, 100, 200, 400, 800, 1600]
for ne in elemRange:
  cntlLines[3] = str(ne) + "             / mesh size (imax)"
  cntlNew = "\n".join(cntlLines)
  cntlOut = open("setflow.cntl","w")
  cntlOut.write(cntlNew)
  cntlOut.close()
  print("\n Running mesh " + str(ne) + "\n")
  subprocess.call(["../../src/multimat1d.exe"])
