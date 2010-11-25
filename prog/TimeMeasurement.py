# -*- coding: utf-8 -*-

#!/usr/bin/env python
import os
import sys
import shutil

res = 2
res2 = 8

#compile

txt1 = """#ifndef _GLOBALSH_
#define _GLOBALSH_

#define NC 3
#define NSPIN 4
#define NDIM 4

//it could be a good idea to define Nt and Ns at compile time
//usually you stick to one volume for quite a while anyways...
#define NSPACE """

txt2 = """
#define VOLSPACE NSPACE*NSPACE*NSPACE
#define VOL4D VOLSPACE*NTIME

#define su2_entries 4

#define PI 	3.14159265358979

#endif
"""
bla = os.getcwd()

#executables must exist!

for i in range (0, res2):
 os.system('make clean')
 os.remove('./globaldefs.h')
 f = open('./globaldefs.h', 'w')
 f.write(txt1)
 f.write(str((i+1)*4) + '\n#define NTIME ' + str((i+1)*4))
 f.write(txt2)
 f.close()
 os.system('make')
 os.rename('hmc', 'hmc_' + str((i+1)*4))
 os.chdir(bla + '/TimeTest')
 os.remove('hmc_' + str((i+1)*4))
 os.chdir(bla)
 shutil.move('hmc_' + str((i+1)*4), bla + '/TimeTest')
 
#execute
#input files must already exist!

os.chdir(bla + '/TimeTest')
for i in range (0,res):
 for j in range (6, res2):
  os.system('./hmc_' + str((j+1)*4) + ' input' + str(i+1));
