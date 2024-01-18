import sys
import os
def save_file(filename,cmds):
    print(filename)
    print(cmds)
    file = open(filename,"w")
    file.write(cmds)
    file.close()

ngen, samplefreq,diagnfreq= sys.argv[2:]

head = f"execute {sys.argv[1]};\nlset nst=6 rates=invgamma;\n"
mb_script = head
mb_script += f"mcmc ngen={ngen} samplefreq={samplefreq} printfreq={samplefreq} diagnfreq={diagnfreq};\n"
mb_script += "q;"
filename = "mb_file"
save_file(filename,mb_script)