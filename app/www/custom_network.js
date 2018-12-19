Shiny.addCustomMessageHandler("jsondata",
  function(message){
    var json_data = message;
    console.log(json_data);

    var w = 1280,
      h = 800,
      rx = w / 2,
      ry = h / 2,
      m0,
      rotate = 0,
      freeze = 0;

    var clickedData = [];

    var splines = [];

    var cluster = d3.layout.cluster()
      .size([360, ry - 120])
      .sort(function(a, b) { return d3.ascending(a.key, b.key); });

    var bundle = d3.layout.bundle();

    var line = d3.svg.line.radial()
      .interpolate("bundle")
      .tension(.85)
      .radius(function(d) { return d.y; })
      .angle(function(d) { return d.x / 180 * Math.PI; });

    // Chrome 15 bug: <http://code.google.com/p/chromium/issues/detail?id=98951>
    var div = d3.select("div.row").insert("div")
      .style("top", "100px")
      .style("left", "450px")
      .style("width", w + "px")
      .style("height", w + "px")
      .style("position", "absolute")
      .style("-webkit-backface-visibility", "hidden");

    var svg = div.append("svg:svg")
      .attr("width", w)
      .attr("height", w)
    .append("svg:g")
      .attr("transform", "translate(" + rx + "," + ry + ")");

    svg.append("svg:path")
      .attr("class", "arc")
      .attr("d", d3.svg.arc().outerRadius(ry - 120).innerRadius(0).startAngle(0).endAngle(2 * Math.PI))
      .on("mousedown", mousedown);

    // testing data for console development
    // var myurl = "https://gist.githubusercontent.com/mbostock/1044242/raw/3ebc0fde3887e288b4a9979dad446eb434c54d08/flare.json"
    // var json_data = await $.getJSON(myurl)

    var nodes = cluster.nodes(root(json_data)),
        links = imports(nodes),
        splines = bundle(links);

    var path = svg.selectAll("path.link")
        .data(links)
      .enter().append("svg:path")
        .attr("class", function(d) { return "link source-" + d.source.key + " target-" + d.target.key; })
        .attr("d", function(d, i) { return line(splines[i]); });

    svg.selectAll("g.node")
        .data(nodes.filter(function(n) { return !n.children; }))
      .enter().append("svg:g")
        .attr("class", "node")
        .attr("id", function(d) { return "node-" + d.key; })
        .attr("transform", function(d) { return "rotate(" + (d.x - 90) + ")translate(" + d.y + ")"; })
      .append("svg:text")
        .attr("dx", function(d) { return d.x < 180 ? 8 : -8; })
        .attr("dy", ".31em")
        .attr("text-anchor", function(d) { return d.x < 180 ? "start" : "end"; })
        .attr("transform", function(d) { return d.x < 180 ? null : "rotate(180)"; })
        .text(function(d) { return d.key; })
        .on("mouseover", mouseover)
        .on("mouseout", mouseout)
        .on("click", mouseclick);

    d3.select("input[type=range]").on("change", function() {
      line.tension(this.value / 100);
      path.attr("d", function(d, i) { return line(splines[i]); });
    });

    d3.select(window)
      .on("mousemove", mousemove)
      .on("mouseup", mouseup);

    function mouse(e) {
    return [e.pageX - rx, e.pageY - ry];
    }

    function mousedown() {
    m0 = mouse(d3.event);
    d3.event.preventDefault();
    }

    function mousemove() {
    if (m0) {
      var m1 = mouse(d3.event),
          dm = Math.atan2(cross(m0, m1), dot(m0, m1)) * 180 / Math.PI;
      div.style("-webkit-transform", "translateY(" + (ry - rx) + "px)rotateZ(" + dm + "deg)translateY(" + (rx - ry) + "px)");
    }
    }

    function mouseup() {
    if (m0) {
      var m1 = mouse(d3.event),
          dm = Math.atan2(cross(m0, m1), dot(m0, m1)) * 180 / Math.PI;

      rotate += dm;
      if (rotate > 360) rotate -= 360;
      else if (rotate < 0) rotate += 360;
      m0 = null;

      div.style("-webkit-transform", null);

      svg
          .attr("transform", "translate(" + rx + "," + ry + ")rotate(" + rotate + ")")
        .selectAll("g.node text")
          .attr("dx", function(d) { return (d.x + rotate) % 360 < 180 ? 8 : -8; })
          .attr("text-anchor", function(d) { return (d.x + rotate) % 360 < 180 ? "start" : "end"; })
          .attr("transform", function(d) { return (d.x + rotate) % 360 < 180 ? null : "rotate(180)"; });
    }
    }

    function mouseover(d) {
        svg.selectAll("path.link.target-" + d.key)
            .classed("target", true)
            .each(updateNodes("source", true));

        svg.selectAll("path.link.source-" + d.key)
            .classed("source", true)
            .each(updateNodes("target", true));
    }

    function mouseout(d) {
      if(!freeze & !clickedData.includes(d)){
        svg.selectAll("path.link.source-" + d.key)
            .classed("source", false)
            .each(updateNodes("target", false));

        svg.selectAll("path.link.target-" + d.key)
            .classed("target", false)
            .each(updateNodes("source", false));
      }
      freeze=0;
    }

    function mouseclick(d){
      //freeze = (freeze) ? 0 : 1
      freeze = 1
      if(!clickedData.includes(d)){
        clickedData.push(d)
      }
      else{
        removeClicked(d)
        svg.selectAll("path.link.source-" + d.key)
            .classed("source", false)
            .each(updateNodes("target", false));

        svg.selectAll("path.link.target-" + d.key)
            .classed("target", false)
            .each(updateNodes("source", false));
      }
    }

    function removeClicked(d, fullName){
      for(var i=0; i < clickedData.length; i++){
        if(d.name === clickedData[i].name){
          clickedData.splice(i-1, i+1)
          break
        }
      }
    }

    function updateNodes(name, value) {
    return function(d) {
      if (value) this.parentNode.appendChild(this);
      svg.select("#node-" + d[name].key).classed(name, value);
    };
    }

    function cross(a, b) {
    return a[0] * b[1] - a[1] * b[0];
    }

    function dot(a, b) {
    return a[0] * b[0] + a[1] * b[1];
    }

    // Lazily construct the package hierarchy from class names.
    function root(json_data) {
      var map = {};

      function find(name, data) {
        var node = map[name], i;
        if (!node) {
          node = map[name] = data || {name: name, children: []};
          if (name.length) {
            node.parent = find(name.substring(0, i = name.lastIndexOf(".")));
            node.parent.children.push(node);
            node.key = name.substring(i + 1);
          }
        }
        return node;
      }

      json_data.forEach(function(d) {
        find(d.name, d);
      });

      return map[""];
    }

    // Return a list of imports for the given array of nodes.
    function imports(nodes) {
      var map = {},
          imports = [];

      // Compute a map from name to node.
      nodes.forEach(function(d) {
        map[d.name] = d;
      });

      // For each import, construct a link from the source to target node.
      nodes.forEach(function(d) {
        if (d.imports) d.imports.forEach(function(i) {
          imports.push({source: map[d.name], target: map[i]});
        });
      });

      return imports;
    }

  });
