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

this_file = "NotInBackground"


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
    strains = [[], []]
    for set_num, set_id in enumerate(['background', 'foreground']):
        for _, _, value, _ in helper.STRAIN_SETS:
            new_strains = form.getvalue(value + set_id)
            if type(new_strains) is list:
                strains[set_num] += new_strains
            elif new_strains is not None:
                strains[set_num].append(new_strains)
    # with open('/playpen/preloaded_data.json') as fp:
    #     data = json.load(fp)
    # with open('/playpen/preloaded_colors.json') as fp:
    #     colors = json.load(fp)
    data, colors = tl.not_in_background(strains[0], strains[1])
    # with open('/playpen/preloaded_data.json', 'w+') as fp:
    # json.dump(data, fp, cls=helper.NumpyEncoder)
    # with open('/playpen/preloaded_colors.json', 'w+') as fp:
    # json.dump(colors, fp, cls=helper.NumpyEncoder)
    plot = bokeh.plotting.figure(y_range=bokeh.models.Range1d(start=tl.offsets[-1] + 10e7, end=0),
                                 tools=['pan', 'wheel_zoom', 'box_zoom', 'save', 'reset', 'resize',
                                        bokeh.models.HoverTool(names=['regions'],
                                                               tooltips=[('Sample', '@sample'),
                                                                         ('Proximal start', '@proximal_start'),
                                                                         ('Proximal size', '@width'),
                                                                         ('Distal start', '@distal_start'),
                                                                         ('Distal size', '@height')])])
    plot_chroms(plot, tl)
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
            sample=combo_data[4]
        )))
        plot.rect('x', 'y', 'active_width', 'active_height', color="#" + hex(color)[2:].zfill(6),
                  source=combo_sources[-1], line_alpha=0, fill_alpha=1.0 / len(strains[1]), name='regions')
    radio_callback = bokeh.models.CustomJS(args={'source' + str(i): source for i, source in enumerate(combo_sources)},
                                           code='var sources = [' + ','.join(
                                               'source' + str(i) for i in xrange(len(combo_sources))) + '];' + '''
            var inputs = document.getElementsByClassName('data-radio');
            var data;
            for (var i = 0; i < inputs.length; i++) {
                data = sources[i].get('data');
                if (inputs[i].checked) {
                    data['active_width'] = data['width'];
                    data['active_height'] = data['height'];
                }
                else {
                    data['active_width'] = 'null';
                    data['active_height'] = 'null';
                }
                sources[i].trigger('change');
            }
            return(false);
            ''')
    radio_button = bokeh.models.widgets.Button(label='', callback=radio_callback)
    bokeh.plotting.show(bokeh.io.vform(plot, radio_button))
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


def radio_buttons(panel):
    panel.table()
    strains = ['Dom', 'Mus', 'Cas']
    panel.tr()
    panel.td()
    panel.td.close()
    for distal in strains:
        panel.td(distal)
        panel.td.close()
    panel.tr.close()
    for i, proximal in enumerate(strains):
        panel.tr()
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