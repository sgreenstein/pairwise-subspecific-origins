import numpy as np
import json
from pairwise_origins import twolocus

with open('../sgreens/pairwise_origins/strain_sets.json') as fp:
    STRAIN_SETS = json.load(fp)


class NumpyEncoder(json.JSONEncoder):
    """ Converts numpy dtypes to the native python equivalent to enable json serialization
    """
    def default(self, o):
        if isinstance(o, np.integer):
            return int(o)
        elif isinstance(o, np.floating):
            return float(o)
        elif isinstance(o, np.ndarray):
            return o.tolist()
        else:
            return super(NumpyEncoder, self).default(o)


def open_control(panel, name):
    panel.div(_class="control-group")
    panel.label(_class="control-label")
    panel.add(name)
    panel.label.close()
    panel.div(_class="controls")


def close_control(panel):
    panel.div.close()  # controls
    panel.div.close()  # control-group


def link_css_and_js(panel):
    panel.script(src="https://ajax.googleapis.com/ajax/libs/jquery/2.1.4/jquery.min.js")
    panel.script.close()
    panel.script(type="text/javascript",
                 src="../sgreens/pairwise_origins/chosen/chosen.jquery.min.js")
    panel.script.close()
    panel.script(type="text/javascript")
    panel.add("""
    function resetSelect() {
        $('.chosen').val('').trigger('chosen:updated');
    }
    """)
    panel.script.close()
    panel.add("""<body onload="resetSelect()">""")
    panel.add("""<body onpageshow = "resetSelect()">""")
    panel.link(rel="stylesheet", type="text/css", href="../sgreens/pairwise_origins/chosen/chosen.min.css")
    panel.link(rel="stylesheet", type="text/css", href="../sgreens/bootstrap/css/bootstrap.min.css")
    panel.link(rel="stylesheet", type="text/css", href="../sgreens/pairwise_origins/all.css")


def strain_set_selector(panel, tl, set_id=''):
    validation_func = "hasSelectedStrains" + set_id
    panel.script('''
    function %s (event) {
        var strainInputs = document.getElementsByClassName("chosen %s");
        var strainInput, j;
        for (var i = 0; i < strainInputs.length; i++) {
            strainInput = strainInputs.item(i);
            for (j = 0; j < strainInput.length; j++) {
                if (strainInput[j].selected) {
                    return true;
                }
            }
        }
        alert("No strains selected");
        return false;
    }
    ''' % (validation_func, set_id), type="text/javascript")
    opened_group_div = False
    for group, text, value, strains in STRAIN_SETS:
        # uncomment the following to add in CC reference
        # if group is not None and not opened_group_div:
        #     opened_group_div = True
        #     panel.h4(group)
        #     panel.div(_class="control-group", style="display: block")
        if group is None:  # skip CC reference for now
            panel.br()
            panel.label(_class="control-label")
            panel.add(text)
            panel.label.close()
            panel.div(_class="controls")
            panel.add("""<select data-placeholder=%s name=%s multiple="multiple" class="chosen %s">""" %
                      (text, value + set_id, set_id))
            for strain in strains:
                if tl.is_available(strain):
                    panel.option(strain, value=strain)
            panel.add('</select>')
            panel.div.close()
    if opened_group_div:
        panel.div.close()
    return validation_func


def select_all_buttons(panel):
    panel.script(type="text/javascript")
    panel.add("""
        $("select").on("chosen:showing_dropdown", function(evnt, params) {
        var chosen = params.chosen,
            $dropdown = $(chosen.dropdown),
            $field = $(chosen.form_field);
        if( !chosen.__customButtonsInitilized ) {
            chosen.__customButtonsInitilized = true;
            var contained = function( el ) {
                var container = document.createElement("div");
                container.appendChild(el);
                return container;
            }
            var width = $dropdown.width();
            var opts = chosen.options || {},
                showBtnsTresshold = opts.disable_select_all_none_buttons_tresshold || 0;
                optionsCount = $field.children().length,
                selectAllText = opts.select_all_text || 'All',
                selectNoneText = opts.uncheck_all_text || 'None';
            if( chosen.is_multiple && optionsCount >= showBtnsTresshold ) {
                var selectAllEl = document.createElement("a"),
                    selectAllElContainer = contained(selectAllEl),
                    selectNoneEl = document.createElement("a"),
                    selectNoneElContainer = contained(selectNoneEl);
                selectAllEl.appendChild( document.createTextNode( selectAllText ) );
                selectNoneEl.appendChild( document.createTextNode( selectNoneText ) );
                $dropdown.prepend("<div class='ui-chosen-spcialbuttons-foot' style='clear:both;border-bottom: 1px solid black;'></div>");
                $dropdown.prepend(selectNoneElContainer);
                $dropdown.prepend(selectAllElContainer);
                var $selectAllEl = $(selectAllEl),
                    $selectAllElContainer = $(selectAllElContainer),
                    $selectNoneEl = $(selectNoneEl),
                    $selectNoneElContainer = $(selectNoneElContainer);
                var reservedSpacePerComp = (width - 25) / 2;
                $selectNoneElContainer.addClass("ui-chosen-selectNoneBtnContainer")
                    .css("float", "right").css("padding", "5px 8px 5px 0px")
                    .css("max-width", reservedSpacePerComp+"px")
                    .css("overflow", "hidden");
                $selectAllElContainer.addClass("ui-chosen-selectAllBtnContainer")
                    .css("float", "left").css("padding", "5px 5px 5px 7px")
                    .css("max-width", reservedSpacePerComp+"px")
                    .css("overflow", "hidden");
                $selectAllEl.on("click", function(e) {
                    e.preventDefault();
                    $field.children().prop('selected', true);
                    $field.trigger('chosen:updated');
                    return false;
                });
                $selectNoneEl.on("click", function(e) {
                    e.preventDefault();
                    $field.children().prop('selected', false);
                    $field.trigger('chosen:updated');
                    return false;
                });
            }
        }
    });
    """)
    panel.script.close()


def visualize_genome(data, tl):
    """ Creates the pairwise genome visualization
    :param data: data to visualize
    :param tl: twolocus instance
    """
    print "content-type: text/html\n"
    print '''
<!DOCTYPE html>
<meta charset="utf-8">
<div style="text-align: center">
<svg class="chart" style="display: inline-block;"></svg>
</div>
<div id="slider"></div>
<button id="zoomout">Zoom out</button>
<script src="https://ajax.googleapis.com/ajax/libs/jquery/2.1.3/jquery.min.js"></script>
<link rel="stylesheet" href="https://ajax.googleapis.com/ajax/libs/jqueryui/1.11.4/themes/smoothness/jquery-ui.css">
<script src="https://ajax.googleapis.com/ajax/libs/jqueryui/1.11.4/jquery-ui.min.js"></script>
<script src="//d3js.org/d3.v3.min.js" charset="utf-8"></script>
<script src="../sgreens/pairwise_origins/d3-zoom-pan-extent.js"></script>
<script type=text/javascript>
var all_data = %s;
var chrom_offsets = %s;
var chrom_sizes = %s;
var chrom_names = %s
</script>
<script src="../sgreens/pairwise_origins/pairwiseGenome.js" charset="utf-8"></script>
        ''' % (data, json.dumps(tl.offsets, cls=NumpyEncoder), json.dumps(tl.sizes, cls=NumpyEncoder),
               json.dumps(twolocus.INT_TO_CHROMO[1:-1]))
