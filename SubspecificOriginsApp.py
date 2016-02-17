import os
import json
import numpy as np
import bokeh.plotting
import bokeh.models

import TwoLocusWebHelper as helper

if __name__ != '__main__':
    from pairwise_origins import twolocus
    import WikiApp
    import markup
    from markup import oneliner as element
    import cgAdmin

this_file = os.path.basename(__file__).replace('App.pyc', '')


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
    panel.h3('Set of Samples')
    has_strains = helper.strain_set_selector(panel, tl)
    panel.script(type="text/javascript")
    panel.add("""$(".chosen").chosen()""")
    panel.script.close()
    helper.select_all_buttons(panel)
    panel.br()
    panel.input(type="hidden", name="target", value="%s.originsVisualization" % this_file)
    panel.input(type="submit", name="submit", value="Submit")
    panel.div.close()  # control group
    panel.form.close()
    panel.div.close()
    panel.script('''$("form").submit(function(event) {return %s(event);});''' % has_strains, type="text/javascript")
    return panel


def originsVisualizationResponse(form):
    # print "content-type: text/json\n"
    coarse_cutoff = 1e7
    panel = markup.page()
    plot_file = 'ss_origins.html'
    bokeh.plotting.output_file(plot_file)
    tl = twolocus.TwoLocus('/csbiodata/public/www.csbio.unc.edu/htdocs/sgreens/pairwise_origins/')
    strains = []
    for _, _, value, _ in helper.STRAIN_SETS:
        new_strains = form.getvalue(value)
        if type(new_strains) is list:
            strains += new_strains
        elif new_strains is not None:
            strains.append(new_strains)
    data, colors = tl.pairwise_frequencies(strains)
    absent_regions = tl.absent_regions(strains)
    plot = bokeh.plotting.figure(y_range=bokeh.models.Range1d(start=tl.offsets[-1] + 10e7, end=0),
                                 height=750, width=750,
                                 background_fill='black',
                                 tools=['box_zoom', 'pan', 'wheel_zoom', 'reset', 'save', 'tap', 'resize',
                                        bokeh.models.HoverTool(
                                            names=["chroms"],
                                            tooltips=[('Proximal', '@proximal'), ('Distal', '@distal')])])
    plot.axis.visible = False
    plot.grid.grid_line_color = None
    for combo_regions, color in zip(data, colors):
        region_widths = np.subtract(combo_regions[1], combo_regions[0])
        region_heights = np.subtract(combo_regions[3], combo_regions[2])
        coarse_indices = np.logical_or(region_widths > coarse_cutoff, region_heights > coarse_cutoff)
        region_widths = region_widths[coarse_indices]
        region_heights = region_heights[coarse_indices]
        x_positions = np.add(np.array(combo_regions[0])[coarse_indices], region_widths/2)
        y_positions = np.add(np.array(combo_regions[2])[coarse_indices], region_heights/2)
        plot.rect(x_positions, y_positions, region_widths, region_heights,
                  color="#" + hex(color)[2:].zfill(6), fill_alpha=1.0 / len(strains), line_alpha=0)
        break
    for combo_regions in absent_regions:
        region_widths = np.subtract(combo_regions[1], combo_regions[0])
        region_heights = np.subtract(combo_regions[3], combo_regions[2])
        coarse_indices = np.logical_or(region_widths > coarse_cutoff, region_heights > coarse_cutoff)
        region_widths = region_widths[coarse_indices]
        region_heights = region_heights[coarse_indices]
        x_positions = np.add(np.array(combo_regions[0])[coarse_indices], region_widths/2)
        y_positions = np.add(np.array(combo_regions[2])[coarse_indices], region_heights/2)
        plot.rect(x_positions, y_positions, region_widths, region_heights, color='white')
        break
    # find coarse data
    # draw chromosomes
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
            chrom_data['proximal'].append(twolocus.INT_TO_CHROMO[i+1])
            chrom_data['distal'].append(twolocus.INT_TO_CHROMO[j+1])
    chrom_data_source = bokeh.models.ColumnDataSource(data=chrom_data)
    plot.rect(x='x', y='y', width='width', height='height', fill_alpha=0, line_color='grey', hover_alpha=0.5,
              name="chroms", source=chrom_data_source)
    text_offsets = tl.offsets.copy()[:-1]
    text_offsets[-3] += 5e7
    plot.text(x=-15e7, y=[offset + 10e7 for offset in text_offsets], text=twolocus.INT_TO_CHROMO[1:],
              text_color='white', text_font_size="10pt")
    plot.text(x=text_offsets, y=tl.offsets[-1] + 10e7, text=twolocus.INT_TO_CHROMO[1:], text_color='white',
              text_font_size="8pt")
    bokeh.plotting.show(plot)
    with open(plot_file) as fp:
        panel.add(fp.read())
    return panel
    helper.visualize_genome(data, tl, len(strains))
    # print json.dumps(data, cls=helper.NumpyEncoder)


if __name__ == '__main__':
    import doctest

    doctest.testmod()
