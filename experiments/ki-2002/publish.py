"""Recompute the equilibrium benchmark and publish code and outputs."""
from pathlib import Path
import os, shutil, subprocess, sys
here = Path(__file__).resolve().parent
out = here/'output'
dest = here.parents[1]/'docs'/here.name
env = dict(os.environ)
env.setdefault('MPLCONFIGDIR', str(out/'.matplotlib'))
subprocess.run([sys.executable, str(here/'reproduce.py'), '--output', str(out),
                '--front-model', str(here.parent/'iwasaki-2012'/'front_model.py')],
               env=env, check=True)
dest.mkdir(parents=True, exist_ok=True)
for name in ['equilibrium.py', 'reproduce.py', 'config.json', 'source-audit.json', 'requirements.txt', 'README.md']:
    shutil.copyfile(here/name, dest/name)
for p in out.iterdir():
    if p.is_file() and p.suffix in ['.json', '.csv', '.png', '.svg']:
        shutil.copyfile(p, dest/p.name)
