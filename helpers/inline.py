import re
from collections import defaultdict
from inspect import getsource
from types import FunctionType


def inline(func):
    return func


def get_indent_level(x):
    return (len(x) - len(x.lstrip())) // 4


def if_else(match, to_inline, new_source, state):
    args = match['args']

    try:
        args = eval(args, to_inline)
        condition, positive, negative = args
        literal = True
    except NameError:
        slitted = args.split(',')
        condition = eval(slitted[0], to_inline)
        positive, negative = slitted[1:]
        literal = False

    if condition:
        target = positive
    else:
        target = negative

    if literal:
        target = repr(target)
    new_source.append('    ' + match['new_name'] + ' = ' + target)


def inline_if(match, to_inline, new_source, state):
    condition = eval(match['condition'].strip(), to_inline)
    state.update({
        'skip_block': state['level'] if not condition else False,
        'skip_since_else': condition,
        'dedent': state['level']
    })


def inline_assignment(match, to_inline, new_source, state):

    old_name = match['name']
    new_name = match['new_name']
    target = to_inline[old_name]

    t = '    ' * state['level']

    if isinstance(target, FunctionType):
        target = target.__name__
    # Function in-lining disabled for now as it gets too complicated
    #    s = getsource(target).split('\n')
    #    new_source.append(t + '# ' + target.__name__)
    #    new_source.append(t + s[0].replace(target.__name__, new_name))
    #    if not s[-1]:   # trim empty line at the end
    #        s.pop(-1)
    #    for l in s[1:]:
    #        new_source.append(t + l)
    s = t + new_name + ' = ' + str(target)
    new_source.append(s)


def inline_return(match, to_inline, new_source, state):
    target = to_inline[match['function']]
    t = '    ' * (state['level'] - 1)
    # TODO: match and replace arguments

    if isinstance(target, FunctionType):
        s = getsource(target).split('\n')

        for l in s[1:]:
            new_source.append(t + l)


def switch_else(match, to_inline, new_source, state):
    if state['skip_since_else']:
        state['skip_block'] = state['level']
        state['skip_since_else'] = False
        state['dedent'] = state['level']
    return False


def compile_with_inline(base_function, name, to_inline, to_pass):

    source = getsource(base_function)
    new_source = []
    lines = source.split('\n')
    state = defaultdict(lambda: False, **{
        'level': 0,
        'global_level': get_indent_level(lines[0]),
    })

    cases = {
        r' *(?P<new_name>.*?) = inline_if_else\((?P<args>.*?)\)$': if_else,
        r' *if inline\((?P<condition>.*?)\):$': inline_if,
        r' *(?P<new_name>.*?) = inline\((?P<name>.*?)\)': inline_assignment,
        r' *else:$': switch_else,
        r' *return inline\((?P<function>.*?)\((?P<arguments>.*?)\)\)$': inline_return
    }

    o = state['global_level'] * 4

    for i, line in enumerate(lines):
        line = line[o:]
        new_level = get_indent_level(line)
        if i == 0 and name:
            line = line.replace(base_function.__name__, name)

        if state['skip_block']:
            if new_level <= state['skip_block']:
                state['skip_block'] = False
                if not line:
                    new_source.append(line)
        else:
            if state['dedent']:
                if new_level <= state['dedent']:
                    state['dedent'] = False
                    state['level'] = state['dedent'] + 1
                line = line[4:]

            for pattern, handler in cases.items():
                match = re.match(pattern, line)
                if match:
                    match = match.groupdict()
                    do_continue = handler(match, to_inline, new_source, state)
                    if not do_continue:
                        break
            else:
                new_source.append(line)

        if line.strip():
            state['level'] = new_level

    source = '\n'.join(new_source)

    f = compile(source, 'module', 'exec')
    module = {**to_pass}
    exec(f, module)
    function = module[name]
    if name:
        function.__name__ = name
    function.__source__ = source
    return function


def inline_if_else(condition, positive, negative):
    return positive if condition else negative
