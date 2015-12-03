import os
import json
import numpy as np
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
    panel.h3('Background Samples')
    has_bg_strains = helper.strain_set_selector(panel, tl, 'background')
    panel.h3('Foreground Samples')
    has_fg_strains = helper.strain_set_selector(panel, tl, 'foreground')
    panel.script("""$(".chosen").chosen()""", type="text/javascript")
    panel.script('''$("form").submit(function () {return %s() && %s();});''' % (has_bg_strains, has_fg_strains),
                 type="text/javascript")
    helper.select_all_buttons(panel)
    panel.br()
    panel.input(type="hidden", name="target", value="%s.uniqueCombos" % this_file)
    panel.input(type="submit", name="submit", value="Submit")
    panel.div.close()  # control group
    panel.form.close()
    panel.div.close()
    return panel


def uniqueCombosResponse(form):
    print "content-type: text/json\n"
    tl = twolocus.TwoLocus('/csbiodata/public/www.csbio.unc.edu/htdocs/sgreens/pairwise_origins/')
    strains = [[], []]
    for set_num, set_id in enumerate(['background', 'foreground']):
        for _, _, value, _ in helper.STRAIN_SETS:
            new_strains = form.getvalue(value + set_id)
            if type(new_strains) is list:
                strains[set_num] += new_strains
            elif new_strains is not None:
                strains[set_num].append(new_strains)
    print '\n'.join(hex(line[0]) + ' ' + ' '.join(map(str, line[1:])) for line in tl.unique_combos(strains[0], strains[1]))
    # print json.dumps(tl.unique_combos(strains[0], strains[1]), cls=helper.NumpyEncoder)
    # panel = markup.page()
    # panel.iframe(src='../sgreens/pairwise_origins/pairwiseGenome.html', width="100%", height="1000px")
    # return panel
    # print "content-type: text/html\n"
    # with open('../sgreens/pairwise_origins/pairwiseGenome.html') as fp:
    #     print fp.read()
    return None


if __name__ == '__main__':
    import doctest
    doctest.testmod()
