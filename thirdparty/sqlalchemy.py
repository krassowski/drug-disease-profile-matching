# The code in this file is based on snippets from documentation of SQLAlchemy:
#     https://docs.sqlalchemy.org/en/latest/orm/extensions/automap.html
# SQLAlchemy and its documentation are licensed under the MIT license:
#     https://github.com/sqlalchemy/sqlalchemy/blob/master/LICENSE
import re
import warnings

import inflect


def name_for_scalar_relationship(base, local_cls, referred_cls, constraint):
    name = referred_cls.__name__.lower()
    local_table = local_cls.__table__
    if name in local_table.columns:
        new_name = name + "_"
        warnings.warn(f'Already detected name {name} present.  using {new_name}')
        return new_name
    return name


_pluralizer = inflect.engine()


def pluralize_collection(base, local_cls, referred_cls, constraint):
    """Produce an 'uncamelized', 'pluralized' class name, e.g.
    'SomeTerm' -> 'some_terms'
    """

    referred_name = referred_cls.__name__
    uncamelized = re.sub(
        r'[A-Z]',
        lambda m: "_%s" % m.group(0).lower(),
        referred_name
    )[1:]
    pluralized = _pluralizer.plural(uncamelized)
    return pluralized


def camelize_classname(base, tablename, table):
    """Produce a 'camelized' class name, e.g.
    'words_and_underscores' -> 'WordsAndUnderscores'
    """

    return str(
        tablename[0].upper()
        +
        re.sub(r'_([a-z])', lambda m: m.group(1).upper(), tablename[1:])
    )


def make_collection_name(base, local_cls, referred_cls, constraint):
    if local_cls is referred_cls:
        return constraint.columns.keys()[0] + "_collection"
