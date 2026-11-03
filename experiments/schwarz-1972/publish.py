"""Recompute and synchronize the Schwarz1972 early-time linear prototype."""
from pathlib import Path
import os
import shutil
import subprocess
import sys
here=Path(__file__).resolve().parent
destination=here.parents[1]/'docs'/'schwarz-1972'
env=dict(os.environ);env.setdefault('MPLCONFIGDIR',str(here/'output'/'.matplotlib'))
subprocess.run([sys.executable,str(here/'reproduce.py')],env=env,check=True)
destination.mkdir(parents=True,exist_ok=True)
for name in ['reproduce.py','config.json','cooling-readings.csv','reference-readings.csv','requirements.txt','source-audit.json','README.md']:
    shutil.copyfile(here/name,destination/name)
for name in ['linear.png','linear.svg','linear.csv','cooling.csv','reference-comparison.csv','validation.json']:
    shutil.copyfile(here/'output'/name,destination/name)
