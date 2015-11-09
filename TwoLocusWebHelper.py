import numpy as np
import json

STRAIN_SETS = [
    ['Classical', 'Classical',
     ["129P1/ReJ", "129P3/J", "129S1SvlmJ", "129S6", "129T2/SvEmsJ", "129X1/SvJ", "A/J", "A/WySnJ", "AEJ/GnLeJ",
      "AEJ/GnRk", "AKR/J", "ALR/LtJ", "ALS/LtJ", "BALB/cByJ", "BALB/cJ", "BDP/J", "BPH/2J", "BPL/1J", "BPN/3J",
      "BTBR T<+>tf/J", "BUB/BnJ", "BXSB/MpJ", "C3H/HeJ", "C3HeB/FeJ", "C57BL/10J", "C57BL/10ScNJ", "C57BL/10ScSnJ",
      "C57BL/6CR", "C57BL/6J", "C57BL/6NCI", "C57BL/6Tc", "C57BLKS/J", "C57BR/cdJ", "C57L/J", "C58/J", "CBA/CaJ",
      "CBA/J", "CE/J", "CHMU/LeJ", "DBA/1J", "DBA/1LacJ", "DBA/2DeJ", "DBA/2HaSmnJ", "DBA/2J", "DDK/Pas",
      "DDY/JclSidSeyFrkJ", "DLS/LeJ", "EL/SuzSeyFrkJ", "FVB/NJ", "HPG/BmJ", "I/LnJ", "IBWSP2", "IBWSR2", "ICOLD2",
      "IHOT1", "IHOT2", "ILS", "ISS", "JE/LeJ", "KK/HlJ", "LG/J", "LP/J", "LT/SvEiJ", "MRL/MpJ", "NOD/ShiLtJ",
      "NON/ShiLtJ", "NONcNZO10/LtJ", "NONcNZO5/LtJ", "NOR/LtJ", "NU/J", "NZB/BlNJ", "NZL/LtJ", "NZM2410/J", "NZO/HlLtJ",
      "NZW/LacJ", "P/J", "PL/J", "PN/nBSwUmabJ", "RF/J", "RHJ/LeJ", "RIIIS/J", "RSV/LeJ", "SB/LeJ", "SEA/GnJ",
      "SEC/1GnLeJ", "SEC/1ReJ", "SH1/LeJ", "SI/Col Tyrp1 Dnahc11/J", "SJL/Bm", "SJL/J", "SM/J", "SSL/LeJ", "ST/bJ",
      "STX/Le", "SWR/J", "TALLYHO/JngJ", "TKDU/DnJ", "TSJ/LeJ", "YBR/EiJ", "ZRDCT Rax<+>ChUmdJ"]],
    ['Wild-derived', 'Wildderived',
     ['22MO', 'BIK/g', 'BULS', 'BUSNA', 'BZO', 'CALB/RkJ', 'CASA/RkJ', 'CAST/EiJ', 'CIM', 'CKN', 'CKS', 'CZECHI/EiJ',
      'CZECHII/EiJ', 'DCA', 'DCP', 'DDO', 'DEB', 'DGA', 'DIK', 'DJO', 'DKN', 'DMZ', 'DOT', 'IS/CamRkJ', 'JF1/Ms',
      'LEWES/EiJ', 'MBK', 'MBS', 'MCZ', 'MDG', 'MDGI', 'MDH', 'MGA', 'MH', 'MOLD/RkJ', 'MOLF/EiJ', 'MOLG/DnJ',
      'MOR/RkJ', 'MPB', 'MSM/Ms', 'PERA/EiJ', 'PERC/EiJ', 'POHN/Deh', 'PWD/PhJ', 'PWK/PhJ', 'RBA/DnJ', 'RBB/DnJ',
      'RBF/DnJ', 'SF/CamEiJ', 'SKIVE/EiJ', 'SOD1/EiJ', 'STLT', 'STRA', 'STRB', 'STUF', 'STUP', 'STUS', 'TIRANO/EiJ',
      'WLA', 'WMP', 'WSB/EiJ', 'ZALENDE/EiJ']],
    ['Wild mice', 'Wild',
     ['BAG102', 'BAG3', 'BAG56', 'BAG68', 'BAG74', 'BAG94', 'BAG99', 'IN13', 'IN17', 'IN25', 'IN34', 'IN38', 'IN40',
      'IN47', 'IN54', 'IN59', 'IN60', 'KCT222', 'MWN1026', 'MWN1030', 'MWN1106', 'MWN1194', 'MWN1198', 'MWN1214',
      'MWN1279', 'MWN1287', 'RDS10105', 'RDS12763', 'RDS13554', 'Yu2095m', 'Yu2097m', 'Yu2099f', 'Yu2113m', 'Yu2115m',
      'Yu2116m', 'Yu2120f']],
    ['CC Founders', 'Founders',
     ["129S1SvlmJ", "A/J", "C57BL/6J", "NOD/ShiLtJ", "NZO/HlLtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ"]],
    ['Sanger', 'Sanger',
     ["129S1SvlmJ", "A/J", "AKR/J", "BALB/cJ", "C3H/HeJ", "CBA/J", "DBA/2J", "LP/J", "NOD/ShiLtJ", "NZO/HlLtJ",
      "PWK/PhJ", "CAST/EiJ", "WSB/EiJ"]]
]


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
    panel.link(rel="stylesheet", type="text/css", href="../sgreens/pairwise_origins/chosen/chosen.min.css")
    panel.link(rel="stylesheet", type="text/css", href="../sgreens/bootstrap/css/bootstrap.min.css")
    panel.link(rel="stylesheet", type="text/css", href="../sgreens/pairwise_origins/all.css")


def strain_set_selector(panel, tl, set_id=''):
    for text, value, strains in STRAIN_SETS:
        panel.br()
        panel.label(_class="control-label")
        panel.add(text)
        panel.label.close()
        panel.div(_class="controls")
        panel.add("""<select data-placeholder=%s name=%s multiple="multiple" class="chosen">""" % (text, value + set_id))
        for strain in strains:
            if tl.is_available(strain):
                panel.option(strain, value=strain)
        panel.add('</select>')
        panel.div.close()


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
