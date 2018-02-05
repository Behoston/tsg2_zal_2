from algorithms.olc import olc_suffix, olc_dynamic, olc
from algorithms.scs import scs
from algorithms.de_bruijn import do_assembly as de_bruijn


algorithms = {
    'SCS': scs,
    'DE_BRUIJN': de_bruijn,
    'OLC': olc,
    'OLC_SUFFIX': olc_suffix,
    'OLC_DYNAMIC': olc_dynamic,
}
