def first(iterable):
    for first in iterable:
        break
    return first


class WarningManager:

    def __init__(self):
        self.history = set()

    def warn_once(self, text):
        if text not in self.history:
            # TODO: use warn and provide custom handler
            print(text)
            self.history.add(text)

    def reset(self):
        self.history.clear()


def on_division_by_zero(fill_with):
    def closure(function):
        def decorated(*args, **kwargs):
            try:
                return function(*args, **kwargs)
            except ZeroDivisionError:
                return fill_with
        decorated.original_function = function
        decorated.decorated_by = f'on_division_by_zero(fill_with={fill_with})'
        decorated.__name__ = function.__name__
        return decorated
    return closure
