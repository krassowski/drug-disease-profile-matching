from tqdm import tqdm_notebook


def progress_bar(*args, leave=False, **kwargs):
    return tqdm_notebook(*args, leave=leave, **kwargs)


tdqm = progress_bar


def no_progress_bar(iterator, total):
    return iterator
