from algorithms.olc import olc_suffix, olc_dynamic
from algorithms.scs import scs

algorithms = {
    'SCS': scs,
    'OLC': olc_suffix,
    'OLC_SUFFIX': olc_suffix,
    'OLC_DYNAMIC': olc_dynamic,
}
