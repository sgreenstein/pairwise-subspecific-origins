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


def parse_position(string):
    """
    >>> parse_position('3,000,000')
    3000000
    >>> parse_position('3M')
    3000000
    >>> parse_position('3kk')
    3000000
    """
    pos = 0
    char_pos = None
    for char_pos, char in enumerate(string):
        if char.isdigit():
            pos *= 10
            pos += int(char)
        elif char != ',':
            break
    # make sure it contained at least 1 digit
    if char_pos is None:
        raise TypeError("Invalid position")
    for char in string[char_pos:]:
        if char.lower() == 'k':
            pos *= 1000
        if char.lower() == 'm':
            pos *= 1000000
    return pos


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
    helper.strain_set_selector(panel, tl)
    for pos_num in ('1', '2'):
        panel.h3('Position ' + pos_num)
        helper.open_control(panel, 'Chromosome')
        panel.select(_class="medium", name="chrom" + pos_num)
        for chrom_num in xrange(1, 19):
            panel.option(chrom_num)
        panel.option('X', value='X')
        panel.select.close()
        helper.close_control(panel)
        helper.open_control(panel, 'Position')
        panel.input(_class="input-medium", name="pos" + pos_num, type="text")
        panel.add("""<p class="help-block">3000000 = 3,000,000 = 3m = 3M = 3000k = 3KK</p>""")
        helper.close_control(panel)
    panel.script(type="text/javascript")
    panel.add("""$(".chosen").chosen()""")
    panel.script.close()
    helper.select_all_buttons(panel)
    panel.br()
    panel.input(type="hidden", name="target", value="%s.countMatrix" % this_file)
    panel.input(type="submit", name="submit", value="Submit")
    panel.div.close()  # control group
    panel.form.close()
    panel.div.close()
    return panel


def countMatrixResponse(form):
    """ Given a set of samples and a pair of loci, find the counts of each combo and the interval bounds
    :return: json encoding of combo counts and interval bounds
    """
    print "content-type: text/json\n"
    tl = twolocus.TwoLocus('/csbiodata/public/www.csbio.unc.edu/htdocs/sgreens/pairwise_origins/')
    strains = []
    for _, value, _ in helper.STRAIN_SETS:
        new_strains = form.getvalue(value)
        if type(new_strains) is list:
            strains += new_strains
        elif new_strains is not None:
            strains.append(new_strains)
    positions = []
    chroms = []
    for pos_num in xrange(2):
        # TODO: error checking
        positions.append(parse_position(form.getvalue('pos' + str(pos_num+1))))
        chroms.append(form.getvalue('chrom' + str(pos_num+1)))
    print json.dumps(tl.sources_at_point_pair(chroms[0], positions[0], chroms[1], positions[1], strains),
                     cls=helper.NumpyEncoder)
    return None


if __name__ == '__main__':
    import doctest
    doctest.testmod()
