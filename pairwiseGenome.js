// constant indices into each array element of data
var COLOR = 4;
var PROX_START = 0;
var PROX_END = 1;
var DIST_START = 2;
var DIST_END = 3;

var color_scale = 20;

var COARSE_CUTOFF = 10000000;  // anything with a dimension smaller than this will be filtered

var genome_length = chrom_offsets[chrom_offsets.length-1];

var chart_height = 800;
var chart_width = 800;

var margin = 17;

var chart = d3.select(".chart")
    .attr("width", chart_width)
    .attr("height", chart_height);

var chart_group = chart.append("g");

var x_axis = chart_group.append("g").attr("id", "x-axis");
var y_axis = chart_group.append("g").attr("id", "y-axis");
var rect_group = chart_group.append("g").attr("id", "rect_group");
var chrom_group = chart_group.append("g").attr("id", "chrom_group");
var x_chrom_labels = chart_group.append("g").attr("id", "x_chrom_labels");
var y_chrom_labels = chart_group.append("g").attr("id", "y_chrom_labels");

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

var translate = d3.scale.linear();
translate.domain([0, genome_length]);
translate.range([margin, chart_width]);
var scale = d3.scale.linear();
scale.domain([0, genome_length]);
scale.range([0, chart_width - margin]);

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
            chart_group.select("#x_chrom_labels > g:nth-child(" + (i+1) + ") > rect")
                .attr("fill", "lightgrey");
            chart_group.select("#y_chrom_labels > g:nth-child(" + (j+1) + ") > rect")
                .attr("fill", "lightgrey");
        })
        .on('mouseout', function (j, i) {
            d3.select(this).attr("stroke", null);
            chart_group.select("#x_chrom_labels > g:nth-child(" + (i+1) + ") > rect")
                .attr("fill", "None");
            chart_group.select("#y_chrom_labels > g:nth-child(" + (j+1) + ") > rect")
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

function drawChromLabels() {
    var groups = x_chrom_labels.selectAll("text")
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
    groups = y_chrom_labels.selectAll("text")
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

function drawAxes(i, j, zoom_scale) {
    var x_ticks = [];
    for (var genome_pos = 0; genome_pos <= chrom_sizes[i]; genome_pos += 300000000 / zoom_scale) {
        x_ticks.push(genome_pos);
    }
    var x_tick_groups = d3.select("#x-axis").selectAll("g").data(x_ticks);
    var new_groups = x_tick_groups.enter().append("g");
    new_groups.append("text").attr("text-anchor", "middle");
    new_groups.append("line");
    x_tick_groups
        .attr("transform", function (d) {
            return "translate(" + translate(d + chrom_offsets[i]) + "," + translate(chrom_offsets[j]) + ")";
        })
        .selectAll("text")
        .text(function (d) {
            return (d / 1000000) + "m";
        })
        .attr("font-size", (16 / zoom_scale).toString());
    x_tick_groups
        .selectAll("line")
        .attr("x1", 0)
        .attr("x2", 0)
        .attr("y1", 0)
        .attr("y2", 1)
        .attr("style", "stroke:#000;stroke-width:0.05");
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
    var zoom = d3.behavior.zoom()
        //.xExtent([x, x + chrom_sizes[i]])// TODO: fix xExtent, yExtent
        //.yExtent([y, y + chrom_sizes[j]])
        .scaleExtent([zoom_scale, 10000])
        .translate([x, y])
        .scale(zoom_scale);
    chart.call(zoom.on("zoom", function () {
        chart_group.attr("transform", "translate(" + d3.event.translate + ")" + "scale(" + d3.event.scale + ")")
    }));
    drawAxes(i, j, zoom_scale);
}

drawRects(coarse_data);
var chromo_rect = drawChroms();
//colorChroms();
drawChromLabels();

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
