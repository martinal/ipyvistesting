from ._version import version_info, __version__

from .volren import *

def _jupyter_nbextension_paths():
    return [{
        'section': 'notebook',
        'src': 'static',
        'dest': 'juypter-volren-widget',
        'require': 'juypter-volren-widget/extension'
    }]
