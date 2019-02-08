from IPython.core.display import clear_output, display

from .namespace import NeatNamespace, HorizontalNamespace, Namespace
from .progress_bar import progress_bar, no_progress_bar
from .table import display_table, bordered_table


def display_one_at_a_time(text):
    clear_output(wait=True)
    display(text)


try:
    from ipywidgets import Output

    class OutputNAtATime:

        def __init__(self, n=5):
            self.n = n
            self.output = Output()
            self.buffer = []

        def activate(self):
            display(self.output)

        def __call__(self, text):
            if len(self.buffer) == self.n:
                self.buffer = self.buffer[1:]
            self.buffer.append(text)
            with self.output:
                clear_output(wait=True)
                display(*self.buffer)

except ImportError:

    class OutputOneAtATime:

        def __call__(self, text):
            return display_one_at_a_time(text)
