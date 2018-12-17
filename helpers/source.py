from functools import reduce

from helpers.gui.source import get_source


def source_for_table(obj, trim_first_decorator=False, type_hints=True):
    import re

    preprocessors = []
    if trim_first_decorator:
        preprocessors.append(lambda source: re.sub('^@(.*?)\n', '', source))
    if not type_hints:
        preprocessors.append(lambda source: re.sub(':(.+?)([=)])', '\\2', source))

    source = get_source(
        obj,
        preprocess=lambda initial_source: (
            reduce(
                lambda source, preprocessor: preprocessor(source),
                [initial_source, *preprocessors]
            )
        )
    )
    source = source._repr_html_().replace('\n', '<br>')
    return f'<span style="text-align: left!important">{source}</span>'
