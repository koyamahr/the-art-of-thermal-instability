"""Recompute the exhibit and synchronize its downloadable files into docs/."""
import os
from pathlib import Path
import shutil
import subprocess
import sys
import argparse

here = Path(__file__).resolve().parent
destination = here.parents[1] / "docs" / "field-1965"
environment = dict(os.environ)
environment.setdefault("MPLCONFIGDIR", str(here / "output" / ".matplotlib"))
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--figure", choices=["1", "2", "3", "all"], default="all")
args = parser.parse_args()
if args.figure in ("1", "all"):
    subprocess.run([sys.executable, str(here / "reproduce-map.py"), "--output", str(here / "output")],
                   env=environment, check=True)
if args.figure in ("2", "all"):
    subprocess.run([sys.executable, str(here / "reproduce.py"), "--output", str(here / "output")],
                   env=environment, check=True)
if args.figure in ("3", "all"):
    subprocess.run([sys.executable, str(here / "reproduce-wave.py"), "--output", str(here / "output")],
                   env=environment, check=True)
destination.mkdir(parents=True, exist_ok=True)
for name in ["reproduce.py", "reproduce-wave.py", "reproduce-map.py", "config.json", "config-wave.json", "config-map.json", "requirements.txt", "environment.txt", "README.md", "README-wave.md", "README-map.md"]:
    shutil.copyfile(here / name, destination / name)
names = []
if args.figure in ("1", "all"):
    names += ["figure-1.svg", "figure-1.png", "figure-1-grid.csv", "figure-1-paths.csv", "validation-map.json"]
if args.figure in ("2", "all"):
    names += ["figure-2.svg", "figure-2.png", "figure-2.csv", "validation.json"]
if args.figure in ("3", "all"):
    names += ["figure-3.svg", "figure-3.png", "figure-3.csv", "validation-wave.json"]
for name in names:
    shutil.copyfile(here / "output" / name, destination / name)
print("Exhibit assets synchronized to", destination)
