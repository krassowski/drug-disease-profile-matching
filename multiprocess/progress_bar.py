"""
MIT licensed, created by Michał Krassowski,
originally published at https://github.com/kn-bibs/pathways-analysis/
"""
from contextlib import contextmanager

from multiprocessing import Queue, Process
from tqdm import tqdm

from .signals import STOP


def progress_bar_worker(queue, total):
    bar = tqdm(total=total)

    for step in iter(queue.get, STOP):
        if step is STOP:
            return
        bar.update(step)


@contextmanager
def progress_bar(total, show=True):
    progress_queue = Queue()

    if show:
        progress = Process(target=progress_bar_worker, args=(progress_queue, total))

        progress.start()

    yield progress_queue

    if show:
        progress_queue.put(STOP)
        progress.join()
