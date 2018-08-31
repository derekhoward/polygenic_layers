import platform
import os
from pathlib import Path

if platform.node() == 'RES-C02TQ1V6G8WL.local':
	RAW_DATA = Path("/Users/derek_howard/projects/HBAsets/data/raw")
	DATA_DIR = Path("/Users/derek_howard/projects/polygenic_layers/data")
	CORES = 6
