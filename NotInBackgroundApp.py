import os
import json
import numpy as np
import bokeh.plotting
import bokeh.models
import twolocus

import TwoLocusWebHelper as helper

if __name__ != '__main__':
    import WikiApp
    import markup
    from markup import oneliner as element
    import cgAdmin

this_file = "NotInBackgroundDev"


def indexPage(form):
    tl = twolocus.TwoLocus('/csbiodata/public/www.csbio.unc.edu/htdocs/sgreens/pairwise_origins/')
    panel = markup.page()
    helper.link_css_and_js(panel)
    panel.div(style="padding:20px 20px;")
    user, permissionLevel, date, time, elapsed = cgAdmin.getCompgenCookie(form)
    editFlag = (form.getvalue("edit") == "True")
    if permissionLevel >= 80:
        panel.add(WikiApp.editTag("%s" % this_file, not editFlag))
    panel.add(WikiApp.getWikiContent("%s" % this_file, editInPlace=(permissionLevel >= 80) and editFlag,
                                     returnPage="./?run=%s" % this_file))
    panel.br()
    panel.form(_class="form-horizontal", action="", method="POST", enctype="multipart/form-data")
    panel.div(_class="control-group")
    panel.h3('Background Samples')
    has_bg_strains = helper.strain_set_selector(panel, tl, 'background')
    panel.h3('Foreground Samples')
    has_fg_strains = helper.strain_set_selector(panel, tl, 'foreground')
    panel.script("""$(".chosen").chosen()""", type="text/javascript")
    panel.script('''$("form").submit(function () {return %s() && %s();});''' % (has_bg_strains, has_fg_strains),
                 type="text/javascript")
    helper.select_all_buttons(panel)
    panel.br()
    panel.input(type="hidden", name="target", value="%s.visualization" % this_file)
    panel.input(type="submit", name="submit", value="Submit")
    panel.div.close()  # control group
    panel.form.close()
    panel.div.close()
    return panel


def visualizationResponse(form):
    panel = markup.page()
    radio_buttons(panel)
    plot_file = 'not_in_background.html'
    bokeh.plotting.output_file(plot_file)
    tl = twolocus.TwoLocus('/csbiodata/public/www.csbio.unc.edu/htdocs/sgreens/pairwise_origins/')
    panel.script(type="text/javascript")
    panel.add('var offsets = ' + json.dumps(tl.offsets, cls=helper.NumpyEncoder) + ';')
    panel.add('var sizes = ' + json.dumps(tl.sizes, cls=helper.NumpyEncoder) + ';')
    panel.add('var chromoToInt = ' + json.dumps(twolocus.CHROMO_TO_INT, cls=helper.NumpyEncoder) + ';')
    panel.add('''
    function getActiveCombo() {
        var inputs = document.getElementsByClassName('data-radio');
        for (var i = 0; i < inputs.length; i++) {
            if (inputs[i].checked) {
                return i;
            }
        }
    }
    var combos = ['Dom', 'Mus', 'Cas'];
    ''')
    panel.script.close()
    strains = [[], []]
    for set_num, set_id in enumerate(['background', 'foreground']):
        for _, _, value, _ in helper.STRAIN_SETS:
            new_strains = form.getvalue(value + set_id)
            if type(new_strains) is list:
                strains[set_num] += new_strains
            elif new_strains is not None:
                strains[set_num].append(new_strains)
    with open('/playpen/preloaded_data.json') as fp:
        data = json.load(fp)
    with open('/playpen/preloaded_colors.json') as fp:
        colors = json.load(fp)
    # data, colors = tl.not_in_background(strains[0], strains[1])
    # with open('/playpen/preloaded_data.json', 'w+') as fp:
    # json.dump(data, fp, cls=helper.NumpyEncoder)
    # with open('/playpen/preloaded_colors.json', 'w+') as fp:
    # json.dump(colors, fp, cls=helper.NumpyEncoder)
    plot = bokeh.plotting.figure(y_range=bokeh.models.Range1d(start=tl.offsets[-1] + 10e7, end=0),
                                 tools=[bokeh.models.HoverTool(names=['chroms'], tooltips=[('Proximal', '@proximal'),
                                                                                           ('Distal', '@distal')])])
    plot.axis.visible = False
    plot.grid.grid_line_color = None
    chrom_source = plot_chroms(plot, tl)
    combo_sources = []
    for i, (combo_data, color) in enumerate(zip(data, colors)):
        width = np.subtract(combo_data[1], combo_data[0])
        height = np.subtract(combo_data[3], combo_data[2])
        combo_sources.append(bokeh.models.ColumnDataSource(data=dict(
            x=np.add(combo_data[0], width / 2),
            y=np.add(combo_data[2], height / 2),
            active_width=width if i == 0 else 'null',
            active_height=height if i == 0 else 'null',
            width=width,
            height=height,
            proximal_start=combo_data[0],
            distal_start=combo_data[2],
            proximal_end=combo_data[1],
            distal_end=combo_data[3],
            sample=combo_data[4]
        )))
        plot.rect('x', 'y', 'active_width', 'active_height', color="#" + hex(color)[2:].zfill(6),
                  source=combo_sources[-1], line_alpha=0, fill_alpha=1.0 / len(strains[1]), name='regions')
    inset_source = bokeh.models.ColumnDataSource(data=dict(
        x=[],
        y=[],
        width=[],
        height=[],
        proximal_start=[],
        proximal_end=[],
        distal_start=[],
        distal_end=[],
        sample=[]
    ))
    inset = bokeh.plotting.figure(y_range=bokeh.models.DataRange1d(flipped=True),
                                  x_range=bokeh.models.DataRange1d(),
                                  tools=['pan', 'wheel_zoom', 'box_zoom', 'save', 'reset', 'resize',
                                         bokeh.models.HoverTool(names=['regions'],
                                                                tooltips=[('Sample', '@sample'),
                                                                          ('Proximal start', '@proximal_start'),
                                                                          ('Proximal end', '@proximal_end'),
                                                                          ('Distal start', '@distal_start'),
                                                                          ('Distal end', '@distal_end')])])
    inset.xaxis.axis_label="Dom"
    inset.yaxis.axis_label="Dom"
    inset.xaxis[0].formatter = bokeh.models.NumeralTickFormatter(format='0.00a')
    inset.yaxis[0].formatter = bokeh.models.NumeralTickFormatter(format='0.00a')
    source_dict = {'source' + str(i): source for i, source in enumerate(combo_sources)}
    source_dict['plot'] = plot
    source_dict['inset'] = inset
    source_dict['chrom_source'] = chrom_source
    source_dict['inset_source'] = inset_source
    radio_callback = bokeh.models.CustomJS(args=source_dict,
                                           code='var sources = [' + ','.join(
                                               'source' + str(i) for i in xrange(len(combo_sources))) + '];' + '''
            var inputs = document.getElementsByClassName('data-radio');
            var data;
            var selected_combo;
            for (var i = 0; i < inputs.length; i++) {
                data = sources[i].get('data');
                if (inputs[i].checked) {
                    data['active_width'] = data['width'];
                    data['active_height'] = data['height'];
                    selected_combo = i;
                }
                else {
                    data['active_width'] = 'null';
                    data['active_height'] = 'null';
                }
                sources[i].trigger('change');
            }
            var distal_chrom = combos[selected_combo%3];
            var proximal_chrom = combos[Math.floor(selected_combo/3)];
            var label = plot.get('left');
            label[0].set('axis_label', distal_chrom);
            label = inset.get('below');
            label[0].set('axis_label', proximal_chrom);
            label = inset.get('left');
            label[0].set('axis_label', distal_chrom);
            label = inset.get('below');
            label[0].set('axis_label', proximal_chrom);
            // clear inset
            inset.set('title', '');
            chrom_source.get('selected')['1d'] = {indices: []};
            inset_data = inset_source.get('data');
            var fields = ['x', 'proximal_start', 'proximal_end', 'distal_start', 'y', 'distal_end', 'width', 'height', 'sample'];
            for (var j = 0; j < fields.length; j++) {
                inset_data[fields[j]] = [];
            }
            inset_source.trigger('change');
            return(false);
            ''')
    chrom_callback = bokeh.models.CustomJS(args=source_dict, code='var sources = [' + ','.join(
        'source' + str(i) for i in xrange(len(combo_sources))) + '];' + '''
    var source = sources[getActiveCombo()];
    var index = cb_obj.get('selected')['1d'].indices;
    var chroms = cb_obj.get('data');
    var proximal = chromoToInt[chroms.proximal[index]];
    var distal = chromoToInt[chroms.distal[index]];
    var i = 0;
    data = source.get('data');
    inset_data = inset_source.get('data');
    var fields = ['x', 'proximal_start', 'proximal_end', 'distal_start', 'y', 'distal_end', 'width', 'height', 'sample'];
    var field, start_index;
    // clear data
    for (var j = 0; j < fields.length; j++) {
        inset_data[fields[j]] = [];
    }
    if (index.length == 0) {
        inset.set('title', '');
        inset_source.trigger('change');
        return(false);
    }

    // set new data
    offset_distal = function (position) {
                        return position - offsets[proximal - 1];
                };
    offset_proximal = function (position) {
                        return position - offsets[distal - 1];
                };
    while (data.proximal_start[i] < offsets[proximal - 1]) {
        i++;  // get through data before proximal chrom
    }
    while (data.proximal_start[i] < offsets[proximal]) {
        while (data.distal_start[i] < offsets[distal - 1] || data.distal_start[i] > offsets[distal]) {
            i++;  // get through data before/after distal chrom
        }
        start_index = i;
        while (data.proximal_start[i] < offsets[proximal] && data.distal_start[i] < offsets[distal]) {
            i++;  // get through data for this chrom pair
        }
        for (var j = 0; j < 3; j++) {
            field = fields[j];
            inset_data[field] = inset_data[field].concat(data[field].slice(start_index, i).map(offset_distal));
        }
        for (var j = 3; j < 6; j++) {
            field = fields[j];
            inset_data[field] = inset_data[field].concat(data[field].slice(start_index, i).map(offset_proximal));
        }
        for (var j = 6; j < 9; j++) {
            field = fields[j];
            inset_data[field] = inset_data[field].concat(data[field].slice(start_index, i));
        }
    }

    // update labels
    var prox_chrom = 'Chr' + chroms.proximal[index];
    var dist_chrom = 'Chr' + chroms.distal[index];
    inset.set('title', prox_chrom + ' x ' + dist_chrom);
    var label = inset.get('left');
    label[0].set('axis_label', label[0].get('axis_label').slice(0, 4) + ' ' + dist_chrom);
    label = inset.get('below');
    label[0].set('axis_label', label[0].get('axis_label').slice(0, 4) + ' ' + prox_chrom);

    // set bounds
    var max_size = Math.max(sizes[distal-1], sizes[proximal-1]);
    var range = inset.get('x_range');
    range.set('bounds', [0, max_size]);
    range.set('start', 0);
    range.set('end', max_size);
    range = inset.get('y_range');
    range.set('bounds', [0, max_size]);
    range.set('start', max_size);
    range.set('end', 0);
    inset_source.trigger('change');
    ''')
    plot.add_tools(bokeh.models.TapTool(names=['chroms'], callback=chrom_callback))
    inset.rect('x', 'y', 'width', 'height', color="black",
               source=inset_source, line_alpha=0, fill_alpha=1.0 / len(strains[1]), name='regions')
    radio_button = bokeh.models.widgets.Button(label='', callback=radio_callback)
    bokeh.plotting.show(bokeh.io.vform(bokeh.io.hplot(plot, inset), radio_button))
    with open(plot_file) as fp:
        panel.add(fp.read())
    button_callbacks(panel)
    return panel


def button_callbacks(panel):
    panel.script(type="text/javascript")
    panel.add('''
    var inputs = document.getElementsByClassName('data-radio');
    for(i=0;i<inputs.length;i++){
      inputs[i].addEventListener("change", function(){
              document.getElementsByClassName('bk-bs-btn')[0].click()
              }, false);
    }
    function hide_buttons() {
        var bokeh_buttons = document.querySelectorAll('button.bk-bs-btn');
        for (var i = 0; i < bokeh_buttons.length; i++) {
            bokeh_buttons[i].style.display = "none";
        }
    }
    document.getElementsByTagName('body')[0].addEventListener("mouseover", hide_buttons, false);
    ''')
    panel.script.close()


def plot_chroms(plot, tl):
    chrom_data = dict(
        x=[],
        y=[],
        width=[],
        height=[],
        proximal=[],
        distal=[],
    )
    for i in xrange(len(tl.sizes)):
        for j in xrange(i, len(tl.sizes)):
            chrom_data['x'].append(tl.offsets[i] + tl.sizes[i] / 2)
            chrom_data['y'].append(tl.offsets[j] + tl.sizes[j] / 2)
            chrom_data['width'].append(tl.sizes[i])
            chrom_data['height'].append(tl.sizes[j])
            chrom_data['proximal'].append(twolocus.INT_TO_CHROMO[i + 1])
            chrom_data['distal'].append(twolocus.INT_TO_CHROMO[j + 1])
    chrom_data_source = bokeh.models.ColumnDataSource(data=chrom_data)
    plot.rect(x='x', y='y', width='width', height='height', fill_alpha=0, line_color='grey', hover_alpha=0.5,
              name="chroms", source=chrom_data_source)
    text_offsets = tl.offsets.copy()[:-1]
    text_offsets[-3] += 5e7
    plot.text(x=-15e7, y=[offset + 10e7 for offset in text_offsets], text=twolocus.INT_TO_CHROMO[1:],
              text_font_size="10pt")
    plot.text(x=text_offsets, y=tl.offsets[-1] + 10e7, text=twolocus.INT_TO_CHROMO[1:],
              text_font_size="8pt")
    return chrom_data_source


def radio_buttons(panel):
    panel.table()
    strains = ['Dom', 'Mus', 'Cas']
    panel.tr()
    panel.td()
    panel.td.close()
    panel.td()
    panel.td.close()
    panel.td(colspan=3)
    panel.add('Distal')
    panel.td.close()
    panel.tr.close()
    panel.tr()
    panel.td()
    panel.td.close()
    panel.td()
    panel.td.close()
    for distal in strains:
        panel.td(distal)
        panel.td.close()
    panel.tr.close()
    for i, proximal in enumerate(strains):
        panel.tr()
        panel.td()
        if i == 1:
            panel.add('Proximal')
        panel.td.close()
        panel.td(proximal)
        panel.td.close()
        for j, distal in enumerate(strains):
            panel.td()
            if i == j == 0:
                panel.input(type="radio", _class="data-radio", name="data-radio", checked=None)
            else:
                panel.input(type="radio", _class="data-radio", name="data-radio")
            panel.td.close()
        panel.tr.close()
    panel.table.close()


if __name__ == '__main__':
    import doctest

    doctest.testmod()
