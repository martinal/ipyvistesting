import ipywidgets as widgets
from traitlets import Unicode, Integer, Float, Bool, Bytes
import pythreejs


npm_module_name = 'juypter-volren-widget'


# TODO: Initially use a robust version of Array with to_json
# TODO: Later optimize using views properly
Float32Buffer = Bytes
Int32Buffer = Bytes


@widgets.register('volrenwidget.VolRenMaterial')
class VolRenMaterial(pythreejs.Material):
    """A prototype volume rendering widget."""
    _view_name = Unicode('VolRenMaterialView').tag(sync=True)
    _model_name = Unicode('VolRenMaterialModel').tag(sync=True)
    _view_module = Unicode(npm_module_name).tag(sync=True)
    _model_module = Unicode(npm_module_name).tag(sync=True)

    # Shaders
    vertexShader = Unicode().tag(sync=True)
    fragmentShader = Unicode().tag(sync=True)

    # Uniforms
    # TODO: Allow fully customized uniforms dict?
    #uniforms = Dict().tag(sync=True)
    time = Float(0.0).tag(sync=True)
    f_min = Float(0.0).tag(sync=True)
    f_max = Float(1.0).tag(sync=True)


@widgets.register('volrenwidget.VolRenGeometry')
class VolRenGeometry(pythreejs.Geometry):
    """A prototype volume rendering widget."""
    _view_name = Unicode('VolRenGeometryView').tag(sync=True)
    _model_name = Unicode('VolRenGeometryModel').tag(sync=True)
    _view_module = Unicode(npm_module_name).tag(sync=True)
    _model_module = Unicode(npm_module_name).tag(sync=True)

    # TODO: Hold updates to rendering while this is true
    #hold = Bool(True).tag(sync=True)

    # Geometry
    vertices = Float32Buffer().tag(sync=True)
    triangles = Int32Buffer().tag(sync=True)

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
