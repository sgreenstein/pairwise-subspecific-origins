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
    print "content-type: text/json\n"
    tl = twolocus.TwoLocus('/csbiodata/public/www.csbio.unc.edu/htdocs/sgreens/pairwise_origins/')
    strains = []
    for _, _, value, _ in helper.STRAIN_SETS:
        new_strains = form.getvalue(value)
        if type(new_strains) is list:
            strains += new_strains
        elif new_strains is not None:
            strains.append(new_strains)
    print json.dumps(tl.pairwise_frequencies(strains), cls=helper.NumpyEncoder)


if __name__ == '__main__':
    import doctest

    doctest.testmod()
