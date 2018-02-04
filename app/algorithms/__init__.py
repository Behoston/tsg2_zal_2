from algorithms.olc import olc_suffix, olc_dynamic, olc
from algorithms.scs import scs

algorithms = {
    'SCS': scs,
    'OLC': olc,
    'OLC_SUFFIX': olc_suffix,
    'OLC_DYNAMIC': olc_dynamic,
}
