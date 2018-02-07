from algorithms.olc import olc_suffix, olc_dynamic, olc
from algorithms.scs import scs as scs
from algorithms.scs import scs as greedy_scs
from algorithms.de_bruijn import do_assembly as de_bruijn


algorithms = {
    'SCS': scs,
    'GREEDY_SCS': greedy_scs,
    'DE_BRUIJN': de_bruijn,
    'OLC': olc,
    'OLC_SUFFIX': olc_suffix,
    'OLC_DYNAMIC': olc_dynamic,
}
