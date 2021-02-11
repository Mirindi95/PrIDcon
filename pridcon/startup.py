import os
from pathlib import Path
import sys 

sys.setrecursionlimit(3000) 

# Default Directory Paths
home_dir = str(Path.home())
PROJECT_DIR = os.path.join(home_dir, ".project03")
CACHE_PATH = os.path.join(PROJECT_DIR, "cache")
LOG_DIR = os.path.join(PROJECT_DIR, "logs")

os.makedirs(CACHE_PATH, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)
