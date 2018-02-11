from contextlib import contextmanager
from time import time


@contextmanager
def timing():
    start = time()
    yield
    print(f"Time consumed {time() - start}")
