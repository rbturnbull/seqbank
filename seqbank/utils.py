# -*- coding: future_typing -*-

from typing import List, Set
from pathlib import Path
import hashlib

def parse_filter(filter:Path|list|set|None=None):
    if not filter:
        return None
    
    if isinstance(filter, Path):
        return set(filter.read_text().strip().split("\n"))
    
    return set(filter)
    

