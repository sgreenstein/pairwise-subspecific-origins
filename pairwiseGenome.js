// constant indices into each array element of data
var COLOR = 4;
var PROX_START = 0;
var PROX_END = 1;
var DIST_START = 2;
var DIST_END = 3;

var color_scale = 20;

var COARSE_CUTOFF = 10000000;  // anything with a dimension smaller than this will be filtered

var genome_length = chrom_offsets[chrom_offsets.length-1];

var chart_height = 1000;
var chart_width = 1500;

var margin = 17;

var chart = d3.select(".chart")
    .attr("width", chart_width)
    .attr("height", chart_height);

var chart_group = chart.append("g");

var rect_group = chart_group.append("g").attr("id", "rect_group");
var chrom_group = chart_group.append("g").attr("id", "chrom_group");
var x_axis = chart_group.append("g").attr("id", "x-axis");
var y_axis = chart_group.append("g").attr("id", "y-axis");

var visible_data;

if (is_ss_origins) {
    visible_data = all_data[0];
}
else {
    visible_data = all_data;
}

var coarse_data;
function filter_coarse_data () {
    coarse_data = visible_data.filter(function (d) {
        return (d[PROX_END] - d[PROX_START] > COARSE_CUTOFF) && (d[DIST_END] - d[DIST_START] > COARSE_CUTOFF)
    });
}
filter_coarse_data();

function dataForChromPair(prox_target, dist_target) {
    return visible_data.filter(function (d) {
            return (d[PROX_START] >= chrom_offsets[prox_target] && d[PROX_END] <= chrom_offsets[prox_target+1] &&
                d[DIST_START] >= chrom_offsets[dist_target] && d[DIST_END] <= chrom_offsets[dist_target+1]);
            });
}

function unflattenIndex(k) {
    return [Math.floor(k / chrom_sizes.length), k % chrom_sizes.length];
}

var chrom_pair_offsets = [];
// initialize data structure for each chrom pair
// first two dimensions are flattened
for (var i = 0; i < chrom_sizes.length; i++) {
    for (var j = i; j < chrom_sizes.length; j++) {
        chrom_pair_offsets.push([chrom_offsets[i], chrom_offsets[j]]);
    }
}

function hexColorString(num) {
    var color = "00000" + num.toString(16);
    return "#" + color.substr(color.length - 6);
}

function rgb(num) {
    return [num >> 16, (num & 0x00ff00) >> 8, num & 0x0000ff];
}

function rgbColorString(r, g, b) {
    return "rgb(" + r + "," + g + "," + b + ")";
}

function scale(pos) {
    return (pos * (chart_height - margin)) / genome_length;
}

function translate(pos) {
    return (pos * (chart_height - margin)) / genome_length + margin;
}

function drawChroms() {
    return chrom_group.selectAll("g")
        .data(chrom_offsets.slice(0, -1))
        .enter().append("g")
        .selectAll("rect")
        .data(function (d, i) {
            var range = [];
            for (j = 0; j <= i; j++) {
                range.push(i);
            }
            return range;
        })
        .enter().append("rect")
        .attr("transform", function (j, i) {
            return "translate(" + translate(chrom_offsets[i]) + "," + translate(chrom_offsets[j]) + ")";
        })
        .attr("fill-opacity", "0")
        .on('mouseover', function (j, i) {
            d3.select(this).attr("stroke", "black");
            chart_group.select("#x-axis > g:nth-child(" + (i+1) + ") > rect")
                .attr("fill", "lightgrey");
            chart_group.select("#y-axis > g:nth-child(" + (j+1) + ") > rect")
                .attr("fill", "lightgrey");
        })
        .on('mouseout', function (j, i) {
            d3.select(this).attr("stroke", null);
            chart_group.select("#x-axis > g:nth-child(" + (i+1) + ") > rect")
                .attr("fill", "None");
            chart_group.select("#y-axis > g:nth-child(" + (j+1) + ") > rect")
                .attr("fill", "None");
        })
        .on('click', function (j, i) {
            zoomToChromPair(i, j);
        })
        .attr("width", function (j, i) {
            return scale(chrom_sizes[i]);
        })
        .attr("height", function (j, i) {
            return scale(chrom_sizes[j]);
        });
}

function colorChroms() {
    chrom_group.selectAll("rect")
        .attr("fill", function (d, k) {
            var areas = unique_areas[k];
            var chrom_indices = unflattenIndex(k);
            var chrom_area = chrom_sizes[chrom_indices[0]] * chrom_sizes[chrom_indices[1]];
            var max_color = null;
            var max_area = 0;
            var total_unique_area = 0;
            for (var color in areas) {
                if (areas.hasOwnProperty(color)) {
                    var area = areas[color];
                    total_unique_area += area;
                    if (area >= max_area) {
                        max_area = area;
                        max_color = color;
                    }
                }
            }
            var max_rgb = rgb(max_color);
            var intensity = 2 / (1 + Math.exp(-1 * color_scale * total_unique_area / chrom_area)) - 1;
            var whiteness = Math.round((1 - intensity) * 0xff);
            max_rgb = max_rgb.map(function (c) {
                return Math.round(c * intensity)
            });
            return rgbColorString(max_rgb[0] + whiteness, max_rgb[1] + whiteness, max_rgb[2] + whiteness);
        });
}

function drawAxes() {
    var groups = x_axis.selectAll("text")
        .data(chrom_names)
        .enter().append("g")
        .attr("transform", function (d, i) {
            return "translate(" + translate(chrom_offsets[i]) + ",0)";
        });
    groups.append("rect")
        .attr("height", margin)
        .attr("width", function (d, i) {
            return scale(chrom_sizes[i])
        })
        .attr("fill", "None");
    groups.append("text")
        .text(function (d) {return d;})
        .attr("alignment-baseline", "before-edge")
        .attr("text-anchor", "middle")
        .attr("transform", function (d, i) {
            return "translate(" + scale(chrom_sizes[i] / 2) + ",0)";
        });
    groups = y_axis.selectAll("text")
        .data(chrom_names)
        .enter().append("g")
        .attr("transform", function (d, i) {
            return "translate(0," + translate(chrom_offsets[i]) + ")";
        });
    groups.append("rect")
        .attr("width", margin)
        .attr("height", function (d, i) {
            return scale(chrom_sizes[i])
        })
        .attr("fill", "None");
    groups.append("text")
        .text(function (d) {return d;})
        .attr("transform", function (d, i) {
            return "translate(0," + scale(chrom_sizes[i] / 2) + ")";
        });
}

function drawRects(data) {
    var color_function, alpha;
    if (is_ss_origins) {
        color_function = function () {
            return hexColorString(source_colors[$("input:radio:checked").attr("value")]);
        };
        alpha = 1 / num_samples;
    }
    else {
        color_function = function (d) { return hexColorString(d[COLOR])};
        alpha = null;
    }
    var rectangles = rect_group.selectAll("rect")
        .data(data, function (d) { return d[PROX_START].toString() + d[DIST_START]});
    rectangles.enter().append("rect").attr("fill-opacity", alpha);
    rectangles
        .attr("transform", function (d) {
            return "translate(" + translate(d[PROX_START]) + "," + translate(d[DIST_START]) + ")";
        })
        .attr("width", function (d) {
            return scale(d[PROX_END] - d[PROX_START]);
        })
        .attr("height", function (d) {
            return scale(d[DIST_END] - d[DIST_START]);
        })
        .attr("fill", color_function);
    rectangles.exit().remove();
}

function zoomToChromPair(i, j) {
    $("#slider").hide();
    $("#zoomout").show();
    chromo_rect.attr("display", "None");
    drawRects(dataForChromPair(i, j));
    var zoom_scale = Math.min(genome_length / chrom_sizes[i],
        genome_length / chrom_sizes[j]);
    var x = -translate(chrom_offsets[i]) * zoom_scale;
    var y = -translate(chrom_offsets[j]) * zoom_scale;
    chart_group.attr("transform",
        "translate(" + x + "," + y + ")scale(" + zoom_scale + ")");
    console.log(x);
    var zoom = d3.behavior.zoom()
        .xExtent([x, x + chrom_sizes[i]])// TODO: fix xExtent, yExtent
        .yExtent([y, y + chrom_sizes[j]])
        .scaleExtent([zoom_scale, 10000])
        .translate([x, y])
        .scale(zoom_scale);
    chart.call(zoom.on("zoom", function () {
        chart_group.attr("transform", "translate(" + d3.event.translate + ")" + "scale(" + d3.event.scale + ")")
    }));
}

drawRects(coarse_data);
var chromo_rect = drawChroms();
//colorChroms();
drawAxes();

$("#zoomout").hide().click( function () {
        $("#slider").show();
        $("#zoomout").hide();
        chromo_rect.attr("display", null);
        chart_group.attr("transform", null);
        var zoom = d3.behavior.zoom().on('zoom', null);
        chart.call(zoom.on("zoom", null));
        drawRects(coarse_data);
    }
);

$("input:radio").change(
    function (e) {
        var source = e.target.value;
        visible_data = all_data[source];
        filter_coarse_data();
        drawRects(coarse_data)
    }
);

//$("#slider").slider(
//    {
//        value:color_scale,
//        min:1,
//        max:10000,
//        step:10,
//        slide: function(event) {
//            color_scale = $("#slider").slider("value");
//            colorChroms();
//        }
//    }
//);
