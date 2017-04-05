import numpy as np
import ipywidgets as widgets
from traitlets import Unicode, Integer, Float, Bool, Bytes, Any, Dict
from traittypes import Array
import pythreejs
from pythreejs.traits_numpy import shape_constraints, array_serialization


# TODO: Avoid duplicate definition here and in __init__.py:
npm_package_name = 'juypter-volren-widget'
npm_package_version = '0.1.0'


def NdArray(dtype=None, shape=(None,), help=""):
    default_shape = tuple(0 if d is None else d for d in shape)
    default_value = np.zeros(shape=default_shape, dtype=dtype)
    return Array(dtype=dtype,
                 default_value=default_value,
                 help=help
                 ).tag(sync=True, **array_serialization
                 ).valid(shape_constraints(*shape))

if 0:
    # Workaround until memoryview bug is fixed in pythreejs/ipywidgets/wherever it is
    def bytes_array_from_json(value, widget):
        # may need to copy the array if the underlying buffer is readonly
        n = np.frombuffer(value['buffer'], dtype=value['dtype'])
        n.shape = value['shape']
        return n
    def bytes_array_to_json(value, widget):
        return {
            'shape': value.shape,
            'dtype': str(value.dtype),
            'buffer': value.tobytes(order='C')
        }
    bytes_array_serialization = dict(to_json=bytes_array_to_json, from_json=bytes_array_from_json)
    def NdArray(dtype=None, shape=(None,), help=""):
        #default_shape = tuple(0 if d is None else d for d in shape)
        #default_value = np.zeros(shape=default_shape, dtype=dtype)
        return Dict(
            #default_value=default_value,
            #help=help
            ).tag(sync=True, **bytes_array_serialization)


@widgets.register
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


@widgets.register
class VolRenGeometry(pythreejs.Geometry):
    """A prototype volume rendering widget."""
    # Required attributes
    _view_name = Unicode('VolRenGeometryView').tag(sync=True)
    _model_name = Unicode('VolRenGeometryModel').tag(sync=True)
    _view_module = Unicode(npm_package_name).tag(sync=True)
    _model_module = Unicode(npm_package_name).tag(sync=True)
    _view_module_version = Unicode(npm_package_version).tag(sync=True)
    _model_module_version = Unicode(npm_package_version).tag(sync=True)

    # Geometry
    position = NdArray(dtype='float32', shape=(None,3), help="Vertex positions")
    faces = NdArray(dtype='uint32', shape=(None,3), help="Triangle vertex indices")

    # Attributes
    # TODO: Allow fully customized attributes dict?
    #attributes = Dict().tag(sync=True)
    f_front = NdArray(dtype='float32', shape=(None,1), help="Front values of function f.")
    f_back = NdArray(dtype='float32', shape=(None,1), help="Back values of function f.")
    s_front = NdArray(dtype='float32', shape=(None,1), help="Front values of function s.")
    s_back = NdArray(dtype='float32', shape=(None,1), help="Back values of function s.")

    # TODO: Higher level interface for setting things in bulk? Using comms instead?
    def set_foo(self, *args):
        # TODO: Hold and update several at once
        pass
