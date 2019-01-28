Shiny.addCustomMessageHandler("jsondata",
  function(message){
    var json_data = message;

    var df = []

    json_data.forEach(function(d){
      d = {name: d.name[0], imports: d.imports, weights: d.weights,
            datasource: d.datasource, Confidence: d.Confidence[0],
            pathway: d.Pathway[0]};
      df.push(d)
    })

    console.log(df);

    var w = 800,
      h = w,
      rx = w / 2,
      ry = h / 2,
      m0,
      rotate = 0,
      headSpace = 63,
      rectW = w-640,
      rectH = w-640;

    var clickedData = [],
        splines = [],
        childrenData = {},
        childrenArray = [];

    var colorMap = ['#ff0000', '#ff8000', '#00ff80', '#0080ff'];
    var windowFields = ['Gene Name:', ' ', 'Gene Family:', ' ', 'Connections:', 'Confidence:'];

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
    var div = d3.select('#igraphViews').insert("div")
      .attr("class", "d3network")
      .style("top", headSpace + "px")
      //.style("left", "0px")
      .style("width", w + "px")
      .style("height", w + "px")
      .style("position", "absolute")
      .style("-webkit-backface-visibility", "hidden");

    var legend = div.append("svg")
      .attr("width", w)
      .attr("height", w/10)
      .attr("class", "legend");

    var svG = div.append("svg:svg")
      .attr("width", w)
      .attr("height", w)

    var svg = svG.append("svg:g")
      .attr("transform", "translate(" + rx + "," + ry + ")");

    var windowArea = svG.append("svg:g")
      .attr("transform", "translate(640,5)")

    windowArea.append("svg:rect")
      .attr("height", rectH)
      .attr("width", rectW)
      .attr("rx", 5)
      .attr("ry", 5)
      .style("stroke", '#0080ff')
      .style("fill", "none")
      .style("stroke-width", 2);

    var windowText = windowArea.selectAll("text")

    windowText.data(windowFields)
      .enter()
        .append("svg:text")
        .attr("class", "vizText")
        .style("fill", "black")
        .attr("transform", function(d) { return "translate(0,20)"; })
        .attr("dy", function(d,i){return(20*i)})
        .attr("dx", 2)
        .text(function(d) { return d; });

    svg.append("svg:path")
      .attr("class", "arc")
      .attr("d", d3.svg.arc().outerRadius(ry - 120).innerRadius(0).startAngle(0).endAngle(2 * Math.PI))
      .on("mousedown", mousedown);

    // testing data for console development
    // var myurl = "https://gist.githubusercontent.com/mbostock/1044242/raw/3ebc0fde3887e288b4a9979dad446eb434c54d08/flare.json"
    // var json_flare = await $.getJSON(myurl)

    // find all names of parentNames
    function findParents(classes){
      var map = {};

      function find(data){
        map[data.parent.name] = data.parent.name;
      }

      classes.forEach(function(d) {
        find(d);
      })

      return Object.keys(map);
    }

    // Lazily construct the package hierarchy from class names.
    function root(classes) {
      var map = {};

      function find(name, data) {
        var node = map[name], i;
        if (!node) {
          node = map[name] = data || {name: name, children: []};
          if (name.length) {
            node.parent = find(name.substring(0, name.indexOf(".")));
            node.parent.children.push(node);
            node.key = name.substring(name.lastIndexOf(".") + 1);
          }
        }
        return node;
      }

      classes.forEach(function(d) {
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

    var nodes = cluster.nodes(root(df)),
        links = imports(nodes),
        splines = bundle(links),
        parents = findParents(df);

    var path = svg.selectAll("path.link")
        .data(links)
      .enter().append("svg:path")
        .attr("class", function(d) { return "link source-" + d.source.key + " target-" + d.target.key; })
        .attr("d", function(d, i) { return line(splines[i]); });

    var colorMapping = function(d){
      for(i in parents){
        if(d.name.includes(parents[i])){
          d.color = colorMap[i]
          return colorMap[i];
        }
      }
    }

    var allNodes = svg.selectAll("g.node")
        .data(nodes.filter(function(n) { return !n.children; }))
      .enter().append("svg:g")
        .attr("class", "node")
        .attr("id", function(d) { return "node-" + d.key; })
        .attr("parentNode", function(d) { return d.parent.name;})
        .attr("transform", function(d) { return "rotate(" + (d.x - 90) + ")translate(" + d.y + ")"; });


    allNodes.append("svg:text")
        .style("fill", colorMapping)
        .style("font-size", function(){
          return (df.length > 350 ? '7.5px' : '10px')
        })
        .attr("dx", function(d) { return d.x < 180 ? 8 : -8; })
        .attr("dy", ".31em")
        .attr("text-anchor", function(d) { return d.x < 180 ? "start" : "end"; })
        .attr("transform", function(d) { return d.x < 180 ? null : "rotate(180)"; })
        .text(function(d) { return d.key; })
        .on("mouseover", mouseover)
        .on("mouseout", mouseout)
        .on("click", mouseclick);

    svg.on("dblclick", mousedbl);

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
      if(!d.clicked){
        if(clickedData.length === 0 || childrenArray.indexOf(d.name) > -1){
          svg.select("#node-" + d.key)
              .style('font-weight', 'bold');
          setvals(d, true)
        }
        if(clickedData.length === 0){
          d3.selectAll(".vizText").remove()
          geneFamily = d.parent.name
          windowFields = ['Gene Name: ',
                          d.key,
                          'Gene Family: ',
                          geneFamily,
                          'Connections: ' + d.datasource.length,
                          'Confidence: ' + d.Confidence]
          /*if(geneFamily != 'Novel Genes'){
            windowFields = ['Gene Name: ',
                            d.key,
                            'Gene Family: ',
                            geneFamily,
                            'Connections: ' + d.datasource.length,
                            'KEGG Confidence: ' + d.keggConf,
                            'STRING Confidence: ' + d.netConf]
                            //'Average Score: ' + Math.round(d3.mean(d.weights) *100) / 100,
                            //'Sources: ' + Array.from(new Set(d.datasource))];
          }
          else{
            windowFields = ['Gene Name: ',
                            d.key,
                            'Gene Family: ',
                            geneFamily,
                            'Connections: ' + d.datasource.length]
                            //'Average Score: ' + Math.round(d3.mean(d.weights) *100) / 100,
                            //'Sources: ' + Array.from(new Set(d.datasource))];
          }*/
          var textColor = d.color;
          windowText.data(windowFields)
            .enter()
              .append("svg:text")
              .attr("class", "vizText")
              .style("fill", function(d, i){
                if(i === 1 || i === 3){
                  return textColor;
                }
              })
              .attr("transform", function(d) { return "translate(0,20)"; })
              .attr("dy", function(d, i){
                return i*20;
              })
              .attr("dx", function(d, i){
                return (i > 3 || i % 2 === 0) ? 2 : 10;
              })
              .text(function(d){return d;});
        }
        else if(childrenArray.indexOf(d.name) > -1){
          d3.selectAll(".vizText").remove()
          n = clickedData[clickedData.length-1]
          cInd = childrenArray.indexOf(d.name)
          windowFields = ['Linked Gene: ', d.key,
                          'Reference Gene: ', n.key,
                          'Connections: ' + d.weights.length,
                          'Score: ' + n.weights[cInd],
                          'Source: ' + n.datasource[cInd]];
          var linkedColor = d.color;
          var refColor = n.color;

          windowText.data(windowFields)
            .enter()
              .append("svg:text")
              .attr("class", "vizText")
              .style("fill", function(d, i){
                if(i === 1){
                  return linkedColor;
                }
                else if(i ===3){
                  return refColor;
                }
              })
              .attr("transform", function(d) { return "translate(0,20)"; })
              .attr("dy", function(d, i){
                return i*20;
              })
              .attr("dx", function(d, i){
                return (i>3 || i % 2 ===0) ? 2 : 10;
              })
              .text(function(d){return d;});
        }
      }
      else{
          svg.select("#node-" + d.key)
              .style('font-weight', 'bold');
          setvals(d, true)
      }
    }

    function mouseout(d) {
      svg.select("#node-" + d.key)
          .style('font-weight', 'normal');
      var i = clickedData.length>1 ? clickedData.length-1 : 0;
      if(!d.clicked){
        setvals(d, false)
        if(clickedData.length>=1){
          setvals(clickedData[i], true, false)
        }
      }
      else if(clickedData.length === 1){
        setvals(d, true, false)
        setvals(clickedData[i], true, false)
      }
      if(clickedData.length === 0){
        clearVizText()
      }
      else{
        d3.selectAll(".vizText").remove()
        n = clickedData[clickedData.length-1]
        windowFields = ['Gene Name: ',
                        n.key,
                        'Gene Family: ',
                        n.parent.name,
                        'Connections: ' + n.datasource.length,
                        'Confidence: ' + n.Confidence]
                        //'Average Score: ' + Math.round(d3.mean(d.weights) *100) / 100,
                        //'Sources: ' + Array.from(new Set(d.datasource))];
        var textColor = n.color;
        windowText.data(windowFields)
          .enter()
            .append("svg:text")
            .attr("class", "vizText")
            .style("fill", function(d, i){
              if(i === 1 || i === 3){
                return textColor;
              }
            })
            .attr("transform", function(d) { return "translate(0,20)"; })
            .attr("dy", function(d, i){
              return i*20;
            })
            .attr("dx", function(d, i){
              return (i > 3 || i % 2 === 0) ? 2 : 10;
            })
            .text(function(d){return d;});
      }
    }

    function mouseclick(d){
      d.clicked = !d.clicked && true;

      clickedData.push(d)

      if(d.clicked){
        removelinks()

        d3.selectAll(".vizText").remove()
        windowFields = ['Gene Name: ',
                        d.key,
                        'Gene Family: ',
                        d.parent.name,
                        'Connections: ' + d.datasource.length,
                        'Confidence: ' + d.Confidence]
                        //'Average Score: ' + Math.round(d3.mean(d.weights) *100) / 100,
                        //'Sources: ' + Array.from(new Set(d.datasource))];
        var textColor = d.color;
        windowText.data(windowFields)
          .enter()
            .append("svg:text")
            .attr("class", "vizText")
            .style("fill", function(d, i){
              if(i === 1 || i === 3){
                return textColor;
              }
            })
            .attr("transform", function(d) { return "translate(0,20)"; })
            .attr("dy", function(d, i){
              return i*20;
            })
            .attr("dx", function(d, i){
              return (i > 3 || i % 2 === 0) ? 2 : 10;
            })
            .text(function(d){return d;});

        //clickedData = Array.from(new Set(clickedData))

        updateChildren(d)

        childrenArray = childrenData[d.name][0]
      }
      else if(clickedData[clickedData.length-1].name === d.name){
        setvals(d, true, false)
      }
      else{
        clickedData.pop()

        updateChildren(d, true)

        setvals(d, false, false)
      }
    }

    function mousedbl(){
      childrenData = {}
      childrenArray = []
      if(clickedData.length>0){
        for(var i in clickedData){
          setvals(clickedData[i], false, false)
          clickedData[i].clicked = false
          svg.select("#node-" + clickedData[i].key)
              .style('font-weight', 'normal');
        }
        clickedData = []
      }
      clearVizText()
    }

    var colorText = function(d){
      for(i in parents){
        if(d === parents[i]){
          return colorMap[i];
        }
      }
    }

    function clearVizText(){
      d3.selectAll(".vizText").remove()
      windowFields = ['Gene Name:', ' ',
                      'Gene Family:', ' ',
                      'Connections:',
                      'Confidence:']
      windowText.data(windowFields)
        .enter()
          .append("svg:text")
          .attr("class", "vizText")
          .style("fill", "black")
          .attr("transform", function(d) { return "translate(0,20)"; })
          .attr("dy", function(d,i){return(20*i)})
          .attr("dx", 2)
          .text(function(d) { return d; });
    }

    function updateNodes(name, value) {
      return function(d) {
        if (value) this.parentNode.appendChild(this);
        svg.select("#node-" + d[name].key).classed(name, value);
      };
    }

    function updateChildren(d, remove=false){
      childrenArray = []

      if(!remove){
        childrenData = {[d.name]: [d.imports]}
      }
      else{
        delete childrenData[d.name]
      }
    }

    function setvals(d, val, block=true){
      if(block){
        if(svg.selectAll("path.link.target-" + d.key).classed("target") === val){
          return;
        }
      }

      svg.selectAll("path.link.target-" + d.key)
          .classed("target", val)
          .each(updateNodes("source", val));

      svg.selectAll("path.link.source-" + d.key)
          .classed("source", val)
          .each(updateNodes("target", val));
    }

    function removelinks(){
      for(var i = 0; i < clickedData.length-1; i++){
        clickedData[i].clicked = false

        svg.selectAll("path.link.target-" + clickedData[i].key)
          .classed('target', false)
          .each(updateNodes("source", false));

        svg.selectAll("path.link.source-" + clickedData[i].key)
          .classed('source', false)
          .each(updateNodes("target", false));
      }
    }

    function cross(a, b) {
      return a[0] * b[1] - a[1] * b[0];
    }

    function dot(a, b) {
      return a[0] * b[0] + a[1] * b[1];
    }

  });
