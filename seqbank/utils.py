
from typing import List, Set
from pathlib import Path

def parse_filter(filter:Path|list|set|None=None) -> Set[str]:
    """ Returns a set of strings to use as a filter.
     
    Args:   
        filter: A list, set or Path to a file containing a list of strings to use as a filter.
                If it is a Path, it should point to a file containing a list of strings, one per line.
                It can also be None, in which case it returns None.
        
    Returns:
        A set of strings to use as a filter.
    """
    if not filter:
        return None
    
    if isinstance(filter, Path):
        return set(filter.read_text().strip().split("\n"))
    
    return set(filter)
    

