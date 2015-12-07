//var data = [[0, 3013441, 126105206, 154919602, 155064859], [0xff, 3013441, 126105206, 155064859, 155065655]];
//var data = [[0, 0, 5, 2, 10], [0xff, 5, 10, 1, 2]];
//var genome_end = 10;

// constant indices into each array element of data
var COLOR = 0;
var PROX_START = 1;
var PROX_END = 2;
var DIST_START = 3;
var DIST_END = 4;

var color_scale = 20;

var COARSE_CUTOFF = 10000000;  // anything with a dimension smaller than this will be filtered

var coarse_data = data.filter( function (d) {
    return (d[PROX_END] - d[PROX_START] > COARSE_CUTOFF) && (d[DIST_END] - d[DIST_START] > COARSE_CUTOFF)
});

function flattenIndex(i, j) {
    return (i * chrom_sizes.length) + j;
}

function unflattenIndex(k) {
    return [Math.floor(k / chrom_sizes.length), k % chrom_sizes.length];
}

// divide data by chromosome combo
var prox_chrom;
var dist_chrom;
var chromo_pairs = [];
var unique_areas = [];
// initialize data structure for each chrom pair
// first two dimensions are flattened
for (var i = 0; i < Math.pow(chrom_sizes.length - 1, 2); i++) {
    chromo_pairs.push([]);
    unique_areas.push({});
}
prox_chrom = 0;
dist_chrom = 0;
var area = 0;
var d;
for (i = 0; i < data.length; i++) {
    d = data[i];
    while (d[PROX_START] > chrom_offsets[prox_chrom+1]) {
        prox_chrom++;
        dist_chrom = 0;
    }
    while (d[DIST_START] > chrom_offsets[dist_chrom+1]) {
        dist_chrom++;
    }
    chromo_pairs[flattenIndex(prox_chrom, dist_chrom)].push(d);
    area = (d[PROX_END] - d[PROX_START]) * (d[DIST_END] - d[DIST_START]);
    unique_areas[flattenIndex(prox_chrom, dist_chrom)][d[COLOR]] =
        (unique_areas[flattenIndex(prox_chrom, dist_chrom)][d[COLOR]] || 0) + area;
}

function hexColorString(num) {
    var color = "00000" + num.toString(16);
    return "#" + color.substr(color.length - 6);
}

function rgb(num) {
    return [num >> 24, (num & 0x00ff00) >> 12, num & 0x0000ff];
}

function rgbColorString(r, g, b) {
    return "rgb(" + r + "," + g + "," + b + ")";
}

function scale(pos) {
    return (pos * chart_height) / chrom_offsets[chrom_offsets.length-1];
}

function colorChroms() {
    d3.selectAll(".chrom")
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
            var alpha = 2 / (1 + Math.exp(-1 * color_scale * total_unique_area / chrom_area)) - 1;
            var whiteness = Math.round((1 - alpha) * 0x0000ff);
            max_rgb = max_rgb.map(function (c) {
                return Math.round(c * alpha)
            });
            return rgbColorString(max_rgb[0] + whiteness, max_rgb[1] + whiteness, max_rgb[2] + whiteness);
        })
        .attr('myattr', color_scale);
}

var chart_height = 1000;
var chart_width = 1500;

var chart = d3.select(".chart")
    .attr("width", chart_width)
    .attr("height", chart_height)
    .call(d3.behavior.zoom().on("zoom", function () {
        chart.attr("transform", "translate(" + d3.event.translate + ")" + "scale(" + d3.event.scale + ")")
    }))
    .append("g");

var uniquities = chart.selectAll("g")
    .data(coarse_data)
    .enter().append("g")
    .attr("transform", function(d) { return "translate(" + scale(d[PROX_START]) + "," + scale(d[DIST_START]) + ")"; });

uniquities.append("rect")
    .attr("width", function(d) { return scale(d[PROX_END] - d[PROX_START]); })
    .attr("height", function(d) { return scale(d[DIST_END] - d[DIST_START]); })
    .attr("fill", function(d) { return hexColorString(d[COLOR]); });

var chromo_rect = chart.selectAll("g")
    .data(chromo_pairs)
    .enter().append("g")
    .attr("transform", function (d, k) {
        var chrom_indices = unflattenIndex(k);
        return "translate(" + scale(chrom_offsets[chrom_indices[0]]) + "," + scale(chrom_offsets[chrom_indices[1]]) + ")";
    });

chromo_rect.append("rect")
    .attr("width", function (d, k) {
        var chrom_indices = unflattenIndex(k);
        return scale(chrom_sizes[chrom_indices[0]]);
    })
    .attr("height", function (d, k) {
        var chrom_indices = unflattenIndex(k);
        return scale(chrom_sizes[chrom_indices[1]]);
    })
    .attr("class", "chrom");

colorChroms();


$("#slider").slider(
    {
        value:color_scale,
        min:1,
        max:10000,
        step:10,
        slide: function(event) {
            color_scale = $("#slider").slider("value");
            colorChroms();
        }
    }
);
