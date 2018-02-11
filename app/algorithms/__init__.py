from algorithms.olc import olc_suffix, olc_dynamic, olc, olc_naive
from algorithms.scs import scs, greedy_scs
from algorithms.de_bruijn import do_assembly as de_bruijn


algorithms = {
    'SCS': scs,
    'GREEDY_SCS': greedy_scs,
    'DE_BRUIJN': de_bruijn,
    'OLC_NAIVE': olc_naive,
    'OLC': olc,
    'OLC_SUFFIX': olc_suffix,
    'OLC_DYNAMIC': olc_dynamic,
}
