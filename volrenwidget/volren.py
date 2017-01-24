import ipywidgets as widgets
from traitlets import Unicode, Integer, Float, Bytes


npm_module_name = 'juypter-volren-widget'


# TODO: Initially use a robust version of Array with to_json
# TODO: Later optimize using views properly
Float32Buffer = Bytes
Int32Buffer = Bytes


@widgets.register('volrenwidget.VolRenWidget')
class VolRenWidget(widgets.DOMWidget):  # TODO: Use pythreejs and make this a Widget instead
    """A prototype volume rendering widget."""
    _view_name = Unicode('VolRenView').tag(sync=True)
    _model_name = Unicode('VolRenModel').tag(sync=True)
    _view_module = Unicode(npm_module_name).tag(sync=True)
    _model_module = Unicode(npm_module_name).tag(sync=True)

    # No rendering should happen until this is set to True
    ready = Bool(False).tag(sync=True)

    # Shaders
    vertex_shader = Unicode().tag(sync=True)
    fragment_shader = Unicode().tag(sync=True)

    # Geometry
    coordinates = Float32Buffer().tag(sync=True)
    triangles = Int32Buffer().tag(sync=True)

    # Uniforms
    # TODO: Allow fully customized uniforms dict?
    #uniforms = Dict().tag(sync=True)
    f_min = Float(0.0).tag(sync=True)
    f_max = Float(1.0).tag(sync=True)

    # Attributes
    # TODO: Allow fully customized attributes dict?
    #attributes = Dict().tag(sync=True)
    f_front = Float32Buffer().tag(sync=True)
    f_back = Float32Buffer().tag(sync=True)
    s_front = Float32Buffer().tag(sync=True)
    s_back = Float32Buffer().tag(sync=True)

    # TODO: Higher level interface for setting things in bulk
    def set_foo(self, *args):
        # TODO: Hold and update several at once
        pass
