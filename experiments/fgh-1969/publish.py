"""Recompute the FGH/GHF ingredient exhibit and synchronize its public assets."""
from pathlib import Path
import os
import shutil
import subprocess
import sys

here=Path(__file__).resolve().parent
destination=here.parents[1]/'docs'/'fgh-1969'
env=dict(os.environ);env.setdefault('MPLCONFIGDIR',str(here/'output'/'.matplotlib'))
subprocess.run([sys.executable,str(here/'reproduce.py'),'--output',str(here/'output')],env=env,check=True)
subprocess.run([sys.executable,str(here/'cold_branch.py'),'--output',str(here/'output')],env=env,check=True)
subprocess.run([sys.executable,str(here/'two_phase.py')],env=env,check=True)
destination.mkdir(parents=True,exist_ok=True)
for name in ['reproduce.py','config.json','source-audit.json','requirements.txt','README.md','cold_branch.py','config-cold.json','two_phase.py','config-two-phase.json']:
    shutil.copyfile(here/name,destination/name)
for name in ['heating.png','heating.svg','heating.csv','cooling-coefficients.png','cooling-coefficients.svg','cooling-coefficients.csv','validation.json','cold-branch.png','cold-branch.svg','cold-branch.csv','validation-cold.json','two-phase.png','two-phase.svg','two-phase.csv','validation-two-phase.json']:
    shutil.copyfile(here/'output'/name,destination/name)
