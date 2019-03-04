Shiny.addCustomMessageHandler("jsondata1",
  function(message){
    var json_data = message;

    var df = []

    json_data.forEach(function(d){
      d = {name: d.name[0], imports: d.imports, weights: d.weights,
            datasource: d.datasource, Confidence: d.Confidence[0]};
      df.push(d)
    })

    console.log(df);

    var w = 800,
      h = w,
      rx = 380,
      ry = 365,
      m0,
      rotate = 0,
      headSpace = 100,
      rectW = w-640,
      rectH = w-650;

    var clickedData = [],
        splines = [],
        childrenData = {},
        childrenArray = [],
        nClicked = 0;

    var colorMap = ['#ff0000', '#ff8000', '#2ECC40', '#0080ff'];
    var windowFields = ['Gene Name:', ' ', 'Connections:', 'Confidence:'];

    var uniquePathways = function(df){
      var lookup = {};
      var result = [];

      for (var i = 0; i<df.length; i++) {
        var name = df[i].parent.name;

        if(name === 'Novel Genes'){continue;}

        if (!(name in lookup)) {
          lookup[name] = 1;
          result.push(name);
        }
      }
      return(result)
    }

    var cluster = d3.layout.cluster()
      .size([360, ry - 120])
      .sort(function(a, b) { return d3.ascending(a.key, b.key); });

    var bundle = d3.layout.bundle();

    var line = d3.svg.line.radial()
      .interpolate("bundle")
      .tension(.85)
      .radius(function(d) { return d.y; })
      .angle(function(d) { return d.x / 180 * Math.PI; });

    if(d3.select("#graphView1")[0][0].children.length === 1){
      d3.select("#graphView1")[0][0].children[0].remove()
    }

    var div = d3.select("#graphView1")
      .attr("class", "d3network")
      .style("top", headSpace + "px")
      .style("width", w + "px")
      .style("height", w + "px")
      .style("position", "absolute")
      .style("-webkit-backface-visibility", "hidden");

    var svG = div.append("svg:svg")
      .attr("width", w)
      .attr("height", w)

    var legend = svG.append("svg:g")
      .attr("transform", "translate(5,5)")

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

    var dubsMessage = svG.append("svg:g")
      .attr("transform", "translate(640, 650)")

    dubsMessage.append("svg:text")
      .style("font-size", "12px")
      .text('** double click to reset');

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

    // find all names of parentNames
    function findParents(classes){
      var map = {};

      function find(data){
        map[data.parent.name] = data.parent.name;
      }

      classes.forEach(function(d) {
        find(d);
      })

      arr = Object.keys(map);
      arr.push(arr.splice(arr.indexOf("Novel Genes"), 1)[0]);

      return arr;
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
        parents = findParents(df),
        pathwayNames = uniquePathways(df);

    if(pathwayNames.length === 1){
      colorMapt = colorMap.slice(0,1)
      colorMapt.push(colorMap[3])
      colorMap = colorMapt
    }
    else if(pathwayNames.length === 2){
      colorMapt = colorMap.slice(0,2)
      colorMapt.push(colorMap[3])
      colorMap = colorMapt
    }

    var colorMapping = function(d){
      for(i in parents){
        if(d.name.includes(parents[i])){
          d.color = colorMap[i]
          return colorMap[i];
        }
      }
    }

    for(i in df){colorMapping(df[i])}

    legend.selectAll("text")
      .data(pathwayNames)
      .enter()
      .append("svg:text")
      .style("fill", function(d, i) {return colorMap[i];})
      .attr("transform", function(d) {return "translate(0,20)";})
      .attr("dy", function(d,i) {return(20*i)})
      //.attr("dx", 2)
      .text(function(d) { return d; });

    var path = svg.selectAll("path.link")
        .data(links)
      .enter().append("svg:path")
        .attr("class", function(d) { return "link source-" + d.source.key + " target-" + d.target.key; })
        .attr("d", function(d, i) { return line(splines[i]); })
        .style('stroke', function(d) {
          if(d.source.parent.name === 'Novel Genes'){
            return d.target.color;
          }
          else{
            return d.source.color;
          }
        });

    var allNodes = svg.selectAll("g.node")
        .data(nodes.filter(function(n) { return !n.children; }))
      .enter().append("svg:g")
        .attr("class", "node")
        .attr("id", function(d) { return "node-" + d.key; })
        .attr("parentNode", function(d) { return d.parent.name;})
        .attr("transform", function(d) { return "rotate(" + (d.x - 90) + ")translate(" + d.y + ")"; });

    allNodes.append("svg:text")
        .style("fill", function(d){ return d.color;})
        .style("font-size", function(){
          return (Math.min(1/Math.log10(df.length)*15, 14))
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
          dm = Math.atan2(cross(m0, m1), dot(m0, m1)) * 180 / Math.PI,
          rotate1 = rotate ? rotate+dm : dm;

      svg.style("-webkit-transform", "translateX(" + rx + "px)translateY(" + ry + "px)rotateZ(" + rotate1 + "deg)");
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
          setvals(d, true, false, pop=5)
        }
        if(clickedData.length === 0){
          d3.selectAll(".vizText").remove()
          paintWindow(d, 1)
        }
        else if(childrenArray.indexOf(d.name) > -1){
          d3.selectAll(".vizText").remove()
          paintWindow(d, 3);
        }
      }
      else{
          setvals(d, true, pop=5)
          svg.select("#node-" + d.key)
              .style('font-weight', 'bold');
      }
    }

    function mouseout(d) {
      var i = clickedData.length>1 ? clickedData.length-1 : 0;
      if(!d.clicked){
        setvals(d, false);
        svg.select("#node-" + d.key)
            .style('font-weight', 'normal');
        if(clickedData.length>=1){
          setvals(clickedData[i], true, false, pop=5)
        }
      }
      else if(clickedData.length === 1){
        setvals(d, true, false, pop=5)
        setvals(clickedData[i], true, false, pop=5)
      }
      if(clickedData.length === 0){
        clearVizText()
      }
      else{
        d3.selectAll(".vizText").remove()
        paintWindow(d, 2)
      }
    }

    function mouseclick(d){
      d.clicked = !d.clicked && true;

      clickedData.push(d)

      if(d.clicked){
        removelinks()

        setvals(d, true, false, pop=5)

        d3.selectAll(".vizText").remove()
        paintWindow(d, 1);

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

      clickeRs = getClicker();

      Shiny.setInputValue("clickedData", clickeRs);
    }

    function getClicker(){
      nConn = clickedData.length
      if(nConn===0){
        return '{}';
      }
      else if(nConn===1){
        clicker = clickedData[0]
        clickeR = '{"Node1": ["' + clicker.name + '"],' +
                  '"Parent1": ["' + clicker.parent.name + '"],' +
                  '"Name1": ["' + clicker.key + '"],' +
                  '"Connections1": ["' + clicker.imports.length + '"],' +
                  '"Confidence": ["' + clicker.Confidence + '"]}'
        return clickeR;
      }
      else{
        var clickeRs = '['
        for(conn in clickedData){
            clicker = clickedData[conn]
            if(conn < nConn-1){
              console.log(true)
              nextClicker = clickedData[parseInt(conn)+1]
              console.log(nextClicker)
              clickeR = '{"Node1": ["' + clicker.name + '"],' +
                        '"Parent1": ["' + clicker.parent.name + '"],' +
                        '"Name1": ["' + clicker.key + '"],' +
                        '"Connections1": ["' + clicker.imports.length + '"],' +
                        '"Confidence": ["' + clicker.Confidence + '"],' +
                        '"Node2": ["' + nextClicker.name + '"],' +
                        '"Parent2": ["' + nextClicker.parent.name + '"],' +
                        '"Name2": ["' + nextClicker.key + '"],' +
                        '"Connections2": ["' + nextClicker.imports.length + '"],' +
                        '"Weight": ["' + clicker.weights[clicker.imports.indexOf(nextClicker.name)] + '"],' +
                        '"Source": ["' + clicker.datasource[clicker.imports.indexOf(nextClicker.name)] + '"]},'
            }
            else{
              clickeR = '{"Node1": ["' + clicker.name + '"],' +
                        '"Parent1": ["' + clicker.parent.name + '"],' +
                        '"Name1": ["' + clicker.key + '"],' +
                        '"Connections1": ["' + clicker.imports.length + '"],' +
                        '"Confidence": ["' + clicker.Confidence + '"]}]'
            }
            clickeRs = clickeRs + clickeR
          }
        return clickeRs;
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

      Shiny.setInputValue("clickedData", '{}');
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


    function updateChildren(d, remove=false){
      childrenArray = []

      if(!remove){
        childrenData = {[d.name]: [d.imports]}
      }
      else{
        delete childrenData[d.name]
      }
    }

    function setvals(d, val, block=true, pop=0){
      if(block){
        if(svg.selectAll("path.link.target-" + d.key).classed("target") === val){
          return;
        }
      }

      svg.selectAll("path.link.target-" + d.key)
          .classed("target", val)
          .style("stroke-opacity", val*.95 + .05)
          .each(updateNodes("source", val, pop));

      svg.selectAll("path.link.source-" + d.key)
          .classed("source", val)
          .style("stroke-opacity", val*.95 + .05)
          .each(updateNodes("source", val, pop));
    }

    function updateNodes(name, value, pop=0) {
      return function(d) {
        if (value) this.parentNode.appendChild(this);

        lastclicked = clickedData[clickedData.length-1];
        nodeBold = value || lastclicked===d[name] ? 'bold' : 'normal';

        svg.select("#node-" + d[name].key)
          .classed(name, value)
          .style('font-weight', nodeBold)
          .transition(50)
          .attr("transform", function(d) { return "rotate(" + (d.x - 90) + ")translate(" + (d.y + pop*value) + ")"; });
      };
    }

    function removelinks(){
      for(var i = 0; i < clickedData.length-1; i++){
        clickedData[i].clicked = false

        svg.selectAll("path.link.target-" + clickedData[i].key)
          .classed('target', false)
          .style("stroke-opacity", 0.05)
          .each(updateNodes("source", false));

        svg.selectAll("path.link.source-" + clickedData[i].key)
          .classed('source', false)
          .style("stroke-opacity", 0.05)
          .each(updateNodes("target", false,));
      }
    }

    function cross(a, b) {
      return a[0] * b[1] - a[1] * b[0];
    }

    function dot(a, b) {
      return a[0] * b[0] + a[1] * b[1];
    }

    function paintWindow(d, option){
      if(option === 1){
        windowFields = ['Gene Name: ',
                        d.key,
                        'Connections: ' + d.datasource.length,
                        'Confidence: ' + d.Confidence]
        windowFields = [].concat.apply([], windowFields);
        var textColor = d.color
        windowText.data(windowFields)
          .enter()
            .append("svg:text")
            .attr("class", "vizText")
            .style("fill", function(d, i){
              if(i===1){
                return textColor;
              }
            })
            .attr("transform", function(d) { return "translate(0,20)"; })
            .attr("dy", function(d, i){
              return i*20;
            })
            .attr("dx", function(d, i){
              if(i===1){
                return 10;
              }
              else{
                return 2;
              }
            })
            .text(function(d){return d;});
      }
      else if(option === 2){
        n = clickedData[clickedData.length-1]
        windowFields = ['Gene Name: ',
                        n.key,
                        'Connections: ' + n.datasource.length,
                        'Confidence: ' + n.Confidence]
        windowFields = [].concat.apply([], windowFields);
        var textColor = n.color;
        var windowL = windowFields.length;
        windowText.data(windowFields)
          .enter()
            .append("svg:text")
            .attr("class", "vizText")
            .style("fill", function(d, i){
              if(i===1){
                return textColor;
              }
            })
            .attr("transform", function(d) { return "translate(0,20)"; })
            .attr("dy", function(d, i){
              return i*20;
            })
            .attr("dx", function(d, i){
              if(i===1){
                return 10;
              }
              else{
                return 2;
              }
            })
            .text(function(d){return d;});
      }
      else{
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
              else if(i === 3){
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

  });
