"""Regenerate the source-rate audit and synchronize only public artifacts."""
from pathlib import Path
import os
import shutil
import subprocess
import sys
here=Path(__file__).resolve().parent
out=here/'output';destination=here.parents[1]/'docs'/here.name
env=dict(os.environ);env.setdefault('MPLCONFIGDIR',str(out/'.matplotlib'))
subprocess.run([sys.executable,str(here/'reproduce.py'),'--output',str(out)],env=env,check=True)
destination.mkdir(parents=True,exist_ok=True)
for name in ['reproduce.py','config.json','source-audit.json','requirements.txt','README.md']:
    shutil.copyfile(here/name,destination/name)
for name in ['rates-and-cooling.png','rates-and-cooling.svg','table2-differences.png','table2-differences.svg','rotational-rates.csv','cooling.csv','table2-comparison.csv','validation.json']:
    shutil.copyfile(out/name,destination/name)
