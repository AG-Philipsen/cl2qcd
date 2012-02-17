#name of the executable
executable = 'inverter'

#debugging options:
#print more information
debug = 0
#print executable-output to stdout
stdout = 1
#save result file under a different name for backup
backup = 0

input_glob = """#global settings
prec=64
use_gpu=true
num_dev=2

startcondition=cold

#fermion settings
fermact=TWISTEDMASS
kappa=0.05
mu=0.2
corr_dir=3
ThetaT=1.

cgmax=100
startcondition=cold
savefrequency=10
fermact=TWISTEDMASS
use_evenodd=yes

solver=BICGSTAB
"""

#arrays for the different tests, this is not nice, but a quick workaround
#the programm will perform tests with all members of this list if no input-file is given
input_var1 = ["""
#variable settings depending on test
NS=16
""",
"""
#variable settings depending on test
NS=24
""",
"""
#variable settings depending on test
NS=32
""",
"""
#variable settings depending on test
NS=48
"""
]

input_var2 = ["""
#variable settings depending on test
NT=4
""",
"""
#variable settings depending on test
NT=8
""",
"""
#variable settings depending on test
NT=12
""",
"""
#variable settings depending on test
NT=16
""",
"""
#variable settings depending on test
NT=20
""",
"""
#variable settings depending on test
NT=24
""",
"""
#variable settings depending on test
NT=28
""",
"""
#variable settings depending on test
NT=32
"""
]
