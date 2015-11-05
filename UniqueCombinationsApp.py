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
    panel.h3('Samples for Set A')
    helper.strain_set_selector(panel, tl, 'A')
    panel.h3('Samples for Set B')
    helper.strain_set_selector(panel, tl, 'B')
    panel.script(type="text/javascript")
    panel.add("""$(".chosen").chosen()""")
    panel.script.close()
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
    for set_num, set_id in enumerate(['A', 'B']):
        for _, value, _ in helper.STRAIN_SETS:
            new_strains = form.getvalue(value + set_id)
            if type(new_strains) is list:
                strains[set_num] += new_strains
            elif new_strains is not None:
                strains[set_num].append(new_strains)
    print json.dumps(tl.unique_combos(strains[0], strains[1]), cls=helper.NumpyEncoder)
    return None


if __name__ == '__main__':
    import doctest
    doctest.testmod()
