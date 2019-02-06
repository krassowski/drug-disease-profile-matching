from select import select
from types import FunctionType
from typing import Dict


def handle_streams(process, handlers: Dict[str, FunctionType]):
    stream_to_handler = {
        stream: handlers[stream_name]
        for stream, stream_name in [
            (process.stdout, 'out'),
            (process.stderr, 'err')
        ]
        if stream_name in handlers
    }

    descriptor_to_stream = {
        stream.fileno(): stream
        for stream in stream_to_handler
    }

    active_descriptors = list(descriptor_to_stream)

    while True:
        selected_descriptors, _, _ = select(active_descriptors, [], [])

        empty_output = True

        for descriptor in selected_descriptors:
            stream = descriptor_to_stream[descriptor]
            line = stream.readline()
            if line:
                line = line.decode('utf-8')
                handler = stream_to_handler[stream]
                handler(line.rstrip())
                empty_output = False

        if empty_output and process.poll() is not None:
            return
