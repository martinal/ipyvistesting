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
            vertices : undefined,
            triangles : undefined,

            // TODO: Customizable set of attributes
            // Attributes
            f_front : undefined,
            f_back : undefined,
            s_front : undefined,
            s_back : undefined
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
            // TODO: When is this called?
            // TODO: Allow updating everything before recreating material.
            // TODO: Allow updating only specific attribute values without recreating material.

            // Setup triangle index array
            var triangles = new Uint32Array(this.model.get('triangles').buffer);

            // Setup vertex coordinate array
            var vertex_components = 3;  // TODO: Might be 2, configure or detect?
            var vertices = new Float32Array(this.model.get('vertices').buffer);

            // Create geometry and connect data to it
            var geometry = new THREE.BufferGeometry();
            geometry.setIndex(new THREE.BufferAttribute(triangles, 1));
            console.log("  triangles:", triangles);
			geometry.addAttribute('position',
                new THREE.BufferAttribute(vertices, vertex_components));
            console.log("  vertices:", vertices);

            // Try to add some other attributes if defined
            console.log("geometry update: attributes");
            // TODO: Make this more configurable somehow, e.g. {"f_front": ["float32", 1]}
            var attribute_names = ["f_front", "f_back", "s_front", "s_back"];
            for (var name of attribute_names) {
                var dataview = this.model.get(name);
                if (dataview !== undefined) {
                    var arr = new Float32Array(dataview.buffer);
			        geometry.addAttribute(name, new THREE.BufferAttribute(arr, 1));
                }
            }

            // Close up
            console.log("geometry update: finalize");
            geometry.computeBoundingSphere();
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
