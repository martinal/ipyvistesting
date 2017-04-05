define(["jupyter-js-widgets", "underscore", "jupyter-threejs", "ndarray"],
function(widgets, _, p3js, ndarray) {
    // Get the three.js library set by pythreejs
    var THREE = window.THREE;

    // Our package name and version
    var npm_package_name = 'juypter-volren-widget';
    var npm_package_version = '0.1.0';


    // FIXME: GeometryModel class is not exported from pythreejs,
    //        at some point this was a problem but now it seems
    //        to work fine to use WidgetModel instead?
    //var GeometryModel = p3js.GeometryModel;
    var GeometryModel = widgets.WidgetModel;


    /// FIXME: Copied from pythreejs
    var typesToArray = {
        int8: Int8Array,
        int16: Int16Array,
        int32: Int32Array,
        uint8: Uint8Array,
        uint16: Uint16Array,
        uint32: Uint32Array,
        float32: Float32Array,
        float64: Float64Array
    }

    /// FIXME: Copied from pythreejs
    var JSONToArray = function(obj, manager) {
        // obj is {shape: list, dtype: string, array: DataView}
        // return an ndarray object
        return ndarray(new typesToArray[obj.dtype](obj.buffer.buffer), obj.shape);
    }

    /// FIXME: Copied from pythreejs
    var arrayToJSON = function(obj, manager) {
        // serialize to {shape: list, dtype: string, array: buffer}
        return {shape: obj.shape, dtype: obj.dtype, buffer: obj.data}
    }

    /// FIXME: Copied from pythreejs
    var array_serialization = { deserialize: JSONToArray, serialize: arrayToJSON };


    // Our custom material model, this holds shaders and uniforms
    var VolRenMaterialModel = p3js.MaterialModel.extend({
        defaults: _.extend(_.result(this, 'p3js.MaterialModel.prototype.defaults'), {
            // Required Model attributes
            _model_name : 'VolRenMaterialModel',
            _view_name : 'VolRenMaterialView',
            _model_module : npm_package_name,
            _view_module : npm_package_name,
            _model_module_version: npm_package_version,
            _view_module_version: npm_package_version,

            // Additional Model attributes
            // Note: Members here should match model in python code,
            // defaults and expected types are defined there

            // Shaders
            vertexShader: 'void main() {}',
            fragmentShader: 'void main() {}',

            // TODO: Customizable set of uniforms
            // Uniforms
            time  : 0.0,
            f_min : 0.0,
            f_max : 1.0
        })
    });

    // Our custom geometry model, this holds mesh and attributes
    var VolRenGeometryModel = GeometryModel.extend({
        defaults: _.extend(_.result(this, 'GeometryModel.prototype.defaults'), {
            // Required Model attributes
            _model_name : 'VolRenGeometryModel',
            _view_name : 'VolRenGeometryView',
            _model_module : npm_package_name,
            _view_module : npm_package_name,
            _model_module_version: npm_package_version,
            _view_module_version: npm_package_version,

            // Additional Model attributes
            // Note: Members here should match model in python code,
            // defaults and expected types are defined there

            // Geometry
            position: ndarray(new Float32Array(), [0, 3]),
            faces: ndarray(new Uint32Array(), [0, 3]),

            // TODO: Customizable set of attributes
            // Attributes
            f_front : ndarray(new Float32Array(), [0, 1]),
            f_back : ndarray(new Float32Array(), [0, 1]),
            s_front : ndarray(new Float32Array(), [0, 1]),
            s_back : ndarray(new Float32Array(), [0, 1]),
        }),
        serializers: _.extend(GeometryModel.serializers, {
            position: array_serialization,
            faces:    array_serialization,
            f_front:  array_serialization,
            f_back:   array_serialization,
            s_front:  array_serialization,
            s_back:   array_serialization,
        })
    });

    // View object mapping material model to THREE.ShaderMaterial object
    var VolRenMaterialView = p3js.ThreeView.extend({
        update: function() {
            console.log("material update: top");
            // TODO: When is this called?
            // TODO: Allow updating everything before recreating material.
            // TODO: Allow updating only specific uniform values without recreating material.

            // TODO: Deal with textures using p3js.DataTexture model
            /*
            var texture = new THREE.TextureLoader().load( "textures/memorial.png" );
            texture.minFilter = THREE.LinearFilter;
            texture.magFilter = THREE.NearestFilter;
            */

            console.log("material update: uniforms");
            var uniforms = {};
            var uniform_names = ["time", "f_min", "f_max"];
            for (var name of uniform_names) {
                var value = this.model.get(name);
                if (value !== undefined) {
                    uniforms[name] = value;
                }
            }

            console.log("material update: create");
			var material = new THREE.ShaderMaterial({
				uniforms: uniforms,
				vertexShader: this.model.get("vertexShader"),
                fragmentShader: this.model.get("fragmentShader")
            });

            // Close up
            console.log("material update: finalize");
            this.replace_obj(material);
        }
    });

    // View object mapping geometry model to THREE.BufferGeometry object
    var VolRenGeometryView = p3js.ThreeView.extend({
        update: function() {
            // TODO: Implement more efficient partial attribute updates

            var geometry = new THREE.BufferGeometry();

            var position = this.model.get('position');
            console.log("position", position);
            var position_array = new Float32Array(position.buffer.buffer);
            console.log("position_array", position_array);
            geometry.addAttribute('position', new THREE.BufferAttribute(position_array, 3));

            var faces = this.model.get('faces');
            var faces_array = new Uint32Array(faces.buffer.buffer);
            console.log("faces_array", faces_array);
            geometry.setIndex(new THREE.BufferAttribute(faces_array, 1));

            // Try to add some other attributes if defined
            // TODO: Make this more configurable somehow, e.g. {"f_front": ["float32", 1]}
            var attribute_names = []; //"f_front", "f_back", "s_front", "s_back"];
            for (var name of attribute_names) {
                var obj = this.model.get(name);
                if (obj !== undefined) {
                    obj_array = new Float32Array(obj.buffer.buffer);
			        geometry.addAttribute(name, new THREE.BufferAttribute(obj_array, 1));
                }
            }

            //geometry.computeVertexNormals();
            geometry.computeBoundingSphere();

            // TODO: What's this for?
            geometry.attributes.position.needsUpdate = true;
            //geometry.attributes.color.needsUpdate = true;

            this.replace_obj(geometry);
        }
    });

    // Export API
    return {
        VolRenMaterialModel : VolRenMaterialModel,
        VolRenMaterialView : VolRenMaterialView,
        VolRenGeometryModel : VolRenGeometryModel,
        VolRenGeometryView : VolRenGeometryView
    };
});
