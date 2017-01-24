var widgets = require('jupyter-js-widgets');
var _ = require('underscore');
//var THREE = require('three');
var p3js = require(["nbextensions/jupyter-threejs/index"]);


var npm_module_name = 'juypter-volren-widget';


// Custom Model. Custom widgets models must at least provide default values
// for model attributes, including `_model_name`, `_view_name`, `_model_module`
// and `_view_module` when different from the base class.
//
// When serialiazing entire widget state for embedding, only values different from the
// defaults will be specified.
var VolRenModel = widgets.DOMWidgetModel.extend({
    defaults: _.extend({}, widgets.DOMWidgetModel.prototype.defaults, {
        _model_name : 'VolRenModel',
        _view_name : 'VolRenView',
        _model_module : npm_module_name,
        _view_module : npm_module_name,

        // No rendering should happen until this is set to True
        ready : false,

        // Shaders
        vertex_shader : undefined,
        fragment_shader : undefined,

        // Geometry
        coordinates : undefined,
        triangles : undefined,

        // Uniforms
        f_min : 0.0,
        f_max : 1.0,

        // Attributes
        f_front : undefined,
        f_back : undefined,
        s_front : undefined,
        s_back : undefined,
    })
});


// Custom View. Renders the widget model.
var VolRenView = widgets.DOMWidgetView.extend({
    ready : false,

    render: function() {
        this.model.on('change:ready', this.ready_changed, this);
    },

    ready_changed: function() {
        var ready = this.model.get('ready');
        var f_min = this.model.get('f_min');
        var f_max = this.model.get('f_max');
        this.el.textContent = "f in range [" + f_min + ", " + f_max + "], ready = " + ready;
    }
});


module.exports = {
    VolRenModel : VolRenModel,
    VolRenView : VolRenView
};
