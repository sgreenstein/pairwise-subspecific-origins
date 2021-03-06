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
    units = ''
    for index, char in enumerate(string):
        if char.lower() in ['k', 'm']:
            units = string[index:]
            string = string[:index]
            break
    pos = float(string.replace(',', ''))
    for char in units:
        if char.lower() == 'k':
            pos *= 1000
        if char.lower() == 'm':
            pos *= 1000000
    pos = int(pos)
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
    has_strains = helper.strain_set_selector(panel, tl)
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
        panel.add('''<input class="input-medium" name="pos%s" required pattern="\d+[\d,]*\.?\d*[mMkK]*">'''
                  % pos_num)
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
    chromo_sizes = {}
    for string, integer in twolocus.CHROMO_TO_INT.iteritems():
        chromo_sizes[string] = tl.sizes[integer-1]
    panel.script('''
    var chromoSizes = %s;
    ''' % json.dumps(chromo_sizes), type="text/javascript")
    panel.script('''
    function parsePosition(string) {
        string = string.replace(',','').toLowerCase();
        var pos = parseFloat(string);
        var units = '';
        var char
        for (var i = 0; i < string.length; i++) {
            char = string[i];
            if (char == 'm' || char == 'k') {
                pos = parseFloat(string.substring(0, i));
                units = string.substring(i);
                break;
            }
        }
        for (i = 0; i < units.length; i++) {
            char = units[i];
            if (char == 'k') {
                pos *= 1000;
            }
            if (char == 'm') {
                pos *= 1000000;
            }
        }
        return Math.floor(pos);
    }

    function isValidIndex(chromosome, position) {
        return (chromosome in chromoSizes) && (0 <= position) && (position <= chromoSizes[chromosome]);
    }

    $("form").submit(function(event) {
        var chrom, pos;
        isValidInput = true;
        for (var i = 1; i <= 2; i++) {
            chrom = document.getElementsByName('chrom' + i)[0].value;
            pos = parsePosition(document.getElementsByName('pos' + i)[0].value);
            if (!isValidIndex(chrom, pos)) {
                alert("Position " + i + " exceeds chromosome length");
                isValidInput = false;
            }
        }
        return %s(event) && isValidInput;
    });
    ''' % has_strains, type="text/javascript")
    return panel


def countMatrixResponse(form):
    """ Given a set of samples and a pair of loci, find the counts of each combo and the interval bounds
    :return: json encoding of combo counts and interval bounds
    """
    # print "content-type: text/json\n"
    tl = twolocus.TwoLocus('/csbiodata/public/www.csbio.unc.edu/htdocs/sgreens/pairwise_origins/')
    strains = []
    for _, _, value, _ in helper.STRAIN_SETS:
        new_strains = form.getvalue(value)
        if type(new_strains) is list:
            strains += new_strains
        elif new_strains is not None:
            strains.append(new_strains)
    positions = []
    chroms = []
    for pos_num in xrange(2):
        positions.append(parse_position(form.getvalue('pos' + str(pos_num + 1))))
        chroms.append(form.getvalue('chrom' + str(pos_num + 1)))
    # print json.dumps(tl.sources_at_point_pair(chroms[0], positions[0], chroms[1], positions[1], strains),
    #                  cls=helper.NumpyEncoder)
    # print tl.sources_at_point_pair(chroms[0], positions[0], chroms[1], positions[1], strains)
    return drawMatrix(tl.sources_at_point_pair(chroms[0], positions[0], chroms[1], positions[1], strains))


def drawMatrix(data):
    panel = markup.page()
    helper.link_css_and_js(panel)
    full_names = {'dom': 'Domesticus', 'mus': 'Musculus', 'cas': 'Castaneus', 'unk': 'Unknown'}
    panel.table()
    panel.tr()
    panel.td('')
    panel.td("Distal interval: Chromosome {}: {:,} - {:,}".format(*tuple(data['Intervals'][1])))
    panel.tr.close()
    panel.tr()
    panel.td("Proximal interval: Chromosome {}: {:,} - {:,}".format(*tuple(data['Intervals'][0])))
    panel.td()
    panel.table(_class="table table-striped")
    panel.tr()
    panel.th('')
    for subspecies in data['Key']:
        panel.th(full_names[subspecies])
    panel.tr.close()
    for subspecies, samples in zip(data['Key'], data['Samples']):
        panel.tr()
        panel.th(full_names[subspecies])
        for sample_set in samples:
            panel.td(', '.join(sample_set) or '-')
        panel.tr.close()
    panel.table.close()
    panel.td.close()
    panel.tr.close()
    panel.table.close()
    return panel


if __name__ == '__main__':
    import doctest

    doctest.testmod()
