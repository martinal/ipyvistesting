import ipywidgets as widgets
from traitlets import Unicode, Integer, Float, Bool, Bytes, Any
import pythreejs


# TODO: Avoid duplicate definition here and in __init__.py:
npm_package_name = 'juypter-volren-widget'
npm_package_version = '0.1.0'


# FIXME: Replace this with new array support in ipywidgets
# TODO: Initially use a robust version of Array with to_json
# TODO: Later optimize using views properly
Float32Buffer = Bytes
Uint32Buffer = Bytes
#Float32Buffer = Any
#Uint32Buffer = Any


@widgets.register('volrenwidget.VolRenMaterial')
class VolRenMaterial(pythreejs.Material):
    """A prototype volume rendering widget."""
    # Required attributes
    _view_name = Unicode('VolRenMaterialView').tag(sync=True)
    _model_name = Unicode('VolRenMaterialModel').tag(sync=True)
    _view_module = Unicode(npm_package_name).tag(sync=True)
    _model_module = Unicode(npm_package_name).tag(sync=True)
    _view_module_version = Unicode(npm_package_version).tag(sync=True)
    _model_module_version = Unicode(npm_package_version).tag(sync=True)

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
    # Required attributes
    _view_name = Unicode('VolRenGeometryView').tag(sync=True)
    _model_name = Unicode('VolRenGeometryModel').tag(sync=True)
    _view_module = Unicode(npm_package_name).tag(sync=True)
    _model_module = Unicode(npm_package_name).tag(sync=True)
    _view_module_version = Unicode(npm_package_version).tag(sync=True)
    _model_module_version = Unicode(npm_package_version).tag(sync=True)

    # TODO: Hold updates to rendering while this is true
    #hold = Bool(True).tag(sync=True)

    # Geometry
    vertices = Float32Buffer().tag(sync=True)
    triangles = Uint32Buffer().tag(sync=True)

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
