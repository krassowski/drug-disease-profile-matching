from pathlib import Path

from jupyter_helpers.notifications import Notifications

ROOT = Path(__file__).parent.parent
SOUNDS_PATH = ROOT / '.jupyter-sounds'

notifications = Notifications(
    success_audio=f'{SOUNDS_PATH}/beep-07.wav',
    failure_audio=f'{SOUNDS_PATH}/beep-05.wav',
    time_threshold=3,  # seconds
    integration='GNOME'
)
