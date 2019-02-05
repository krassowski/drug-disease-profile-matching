from pathlib import Path
from time import time
from IPython import get_ipython
from IPython.display import Audio, display


SOUNDS_PATH = Path(__file__).parent.parent / '.jupyter-sounds'


class InvisibleAudio(Audio):
    def _repr_html_(self):
        audio = super()._repr_html_()
        audio = audio.replace('<audio', f'<audio onended="this.parentNode.removeChild(this)"')
        return f'<div style="display:none">{audio}</div>'
        
        
def play_sound_on_exception(self, etype, value, tb, tb_offset=None):
    self.showtraceback((etype, value, tb), tb_offset=tb_offset)
    # http://www.soundjay.com/button/beep-05.wav
    display(InvisibleAudio(filename=f'{SOUNDS_PATH}/beep-05.wav', autoplay=True))

    
class Beeper:

    def __init__(self, threshold, **audio_kwargs):
        self.threshold = threshold
        self.start_time = None    # time in sec, or None
        self.audio = audio_kwargs
        
    def pre_execute(self):
        if not self.start_time:
            self.start_time = time()

    def post_execute(self):
        end_time = time()
        if self.start_time and end_time - self.start_time > self.threshold:
            audio = InvisibleAudio(**self.audio, autoplay=True)
            display(audio)
        self.start_time = None

    post_execute.is_beeper_handler = True
    pre_execute.is_beeper_handler = True


def remove_old_beepers(ipython):
    for event, callbacks in ipython.events.callbacks.items():
        for callback in callbacks:
            if hasattr(callback, 'is_beeper_handler'):
                ipython.events.unregister(event, callback)


# http://www.soundjay.com/button/beep-07.wav
beeper = Beeper(1, filename=f'{SOUNDS_PATH}/beep-07.wav')

ipython = get_ipython()
remove_old_beepers(ipython)
ipython.events.register('pre_execute', beeper.pre_execute)
ipython.events.register('post_execute', beeper.post_execute)
ipython.set_custom_exc((Exception,), play_sound_on_exception)
