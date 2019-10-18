Shiny.addCustomMessageHandler("jsondata2",
  function(message){
    var json_data = message;

    var df = []

    json_data.forEach(function(d){
      d = {name: d.name[0], imports: d.imports, weights: d.weights,
            datasource: d.datasource, Confidence: d.Confidence[0], order: d.order[0], color: d.Color[0]};
      df.push(d)
    })

    //console.log(df);

    var w = 800,
      h = w,
      rx = 380,
      ry = 365,
      m0,
      rotate = 0,
      headSpace = 100,
      columnSpace = 250,
      niaidBannerY = 170,
      rectW = w-640,
      rectH = w-650,
      splines = [];

    let clickedData2 = [],
        childrenData2 = {},
        childrenArray2 = [];

    var novelColor = "green",
        color1 = "blue",
        color2 = "red",
        color3 = "saddlebrown",
        color12 = "#a09d01",
        color13 = "#6b5b3e",
        color23 = "#ce6702",
        color123 = "#6d1c8e",
        colorMap = [color1, color2, color3, novelColor, color123, color12, color23, color13],
        windowFields = ['Gene Name:', ' ', 'Interactions:', 'Screen Input:'];

    var cluster = d3.layout.cluster()
      .size([360, ry - 120])
      .sort(function(a, b) { return d3.ascending(a.order, b.order)})

    var bundle = d3.layout.bundle();
    var line = d3.svg.line.radial()
      .interpolate("bundle")
      .tension(.85)
      .radius(function(d) { return d.y; })
      .angle(function(d) { return d.x / 180 * Math.PI; });

    var div = d3.select("#graphView2")
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

    //slider for tension
    div.append('input')
      .attr("type", "range")
      .attr("class", "tension")
      .attr("min", 0)
      .attr("max", 100)
      .on("change", function() {
        line.tension(this.value / 100);
        path.attr("d", function(d, i) { return line(splines[i]); });
      });

    //slider for text size
    div.append('input')
      .attr("type", "range")
      .attr("class", "textSize")
      .attr("min", 1)
      .attr("max", 18)
      .attr("step", 0.2)
      .attr("value", Math.min(1/Math.log10(df.length)*15, 14))
      .on("change", function() {
        d3.selectAll(".node2 text").style("font-size", this.value);
      });

    //slider for node size
    div.append('input')
      .attr("type", "range")
      .attr("class", "nodeSize")
      .attr("min", 0.5)
      .attr("max", 10)
      .attr("step", 0.1)
      .attr("value", Math.min(1/Math.log10(df.length)*5, 4))
      .on("change", function() {
        d3.selectAll(".node2 circle").style("r", this.value)
      });

    //slider for window and legend text size
    div.append('input')
      .attr("type", "range")
      .attr("class", "wTextSize")
      .attr("min", 1)
      .attr("max", 24)
      .attr("step", 1)
      .attr("value", 14)
      .on("change", function() {
        legend.selectAll("text").attr("font-size", this.value)
      });

    sliderText = ['Tension', 'Text Size', 'Node Size', 'Legend Text']

    var sliderSVG = svG.append("svg:g")
      .attr("width", w)
      .attr("height", 20)
      .attr("class", "slidersvg")
      .attr("transform", "translate(700,360)")

    sliderSVG.selectAll("text")
      .data(sliderText)
      .enter()
      .append("svg:text")
      .attr("font-family", "Helvetica")
      .attr("font-size", 14)
      .attr("dy", function(d,i){return i*90})
      .text(function(d){return d;});

    // removes the network div once another network is chosen in shiny
    if(d3.select("#graphView2")[0][0].children.length > 5){
      while(d3.select("#graphView2")[0][0].children.length > 5){
        d3.select("#graphView2")[0][0].children[0].remove()
      }
    }

    // button for going back in clicked pathways
    var revert = svG.append("svg:g")
      .attr("class", "revert")
      .attr("transform", "translate(700,265)")

    revert.append("svg:rect")
      .attr("height", 20)
      .attr("width", 95)
      .attr("rx", 10)
      .attr("ry", 10)
      .style("stroke", "black")
      .style("fill", "none")
      .style("stroke-width", 2);

    revert.append("svg:text")
      .attr("font-family", "Helvetica")
      .style("text-anchor", "center")
      .attr("dy", 15)
      .attr("dx", 8)
      .text('Revert Click')
      .on("click", revertOne);

    // setting up nodes , links, splines, and parents for references and visualizations
    // of network graph
    var nodes = cluster.nodes(root(df)),
        links = imports(nodes),
        splines = bundle(links),
        parents = findParents(df);

    var svg = svG.append("svg:g")
      .attr("class", "networkGraph")
      .attr("transform",  "translate(" + rx + "," + ry + ")");

    var windowArea = svG.append("svg:g")
      .attr("transform", "translate(638,5)")

    windowArea.append("svg:rect")
      .attr("height", rectH)
      .attr("width", rectW)
      .attr("rx", 5)
      .attr("ry", 5)
      .style("stroke", novelColor)
      .style("fill", "none")
      .style("stroke-width", 2);

    var dubsMessage = svG.append("svg:g")
      .attr("transform", "translate(694, 170)")
      .attr("class", "resetButton");

    dubsMessage.append("svg:rect")
      .attr("height", 20)
      .attr("width", 102)
      .attr("rx", 10)
      .attr("ry", 10)
      .style("stroke", "black")
      .style("fill", "none")
      .style("stroke-width", 2);

    dubsMessage.append("svg:text")
      .attr("font-family", "Helvetica")
      .style("text-anchor", "center")
      .attr("dy", 15)
      .attr("dx", 9)
      .text('Click to reset')
      .on("click", mousedbl);

    var windowText = windowArea.selectAll("text")

    windowText.data(windowFields)
      .enter()
        .append("svg:text")
        .attr("class", "vizText")
        .style("fill", "black")
        .attr("transform", function(d) { return "translate(0,20)"; })
        .attr("font-family", "Helvetica")
        .attr("font-size", 14)
        .attr("dy", function(d,i){return(20*i)})
        .attr("dx", 2)
        .text(function(d) { return d; });

    var pathButton = svG.append("svg:g")
      .attr("transform", "translate(715, 200)")
      .attr("class", "pathButton")

    pathButton.append("svg:rect")
      .attr("height", 55)
      .attr("width", 80)
      .attr("rx", 10)
      .attr("ry", 10)
      .style("stroke", "black")
      .style("fill", "none")
      .style("stroke-width", 2);

    pathButtonTxt = ["Highlight", "Clicked", "Pathways"]

    pathButton.selectAll("text")
      .data(pathButtonTxt)
      .enter()
        .append("svg:text")
        .style("text-anchor", "center")
        .attr("dx", 7)
        .attr("dy", function(d, i){return 16*i + 16;})
        .attr("font-family", "Helvetica")
        .text(function(d){return d;})
        .on("click", drawConnections);

    svg.append("svg:path")
      .attr("class", "arc")
      .attr("d", d3.svg.arc().outerRadius(ry - 120).innerRadius(0).startAngle(0).endAngle(2 * Math.PI))

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
      arr.splice(arr.indexOf("Novel"), 1);
      for(var i=0;i<arr.length;i++){
        arr[i]=arr[i] + " (TRIAGE hits)";
      }
      arr.push("Novel")

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

    var colorMapping = function(df){
      newCols = []
      for(i in df){
        if(!(newCols.includes(df[i].color))){
          newCols.push(df[i].color)
        }
      }
      newCols.splice(newCols.indexOf(novelColor), 1);
      newCols.push(novelColor);
      return newCols;
    }

    startingColors = colorMapping(df)

    var ordColors = []

    for(j in colorMap){
      if(startingColors.includes(colorMap[j])){
        ordColors.push(colorMap[j])
      }
    }

    var novelInd = ordColors.indexOf(novelColor)

    var getLegendColors = function(){
      check = 0
      for(i in parents){
        if(parents[i].includes(" & ")){
          check = 1
        }
      }
      if(check){
        lColors = ["#000000"];
        lColors = lColors.concat(ordColors.slice(0, novelInd+1));
        lColors.push("#000000", "#000000", "#000000");
        lColors = lColors.concat(ordColors.slice(novelInd+1));
      }
      else{
        lColors = startingColors;
      }
      return lColors;
    }

    var legendColors = getLegendColors();

    var pathwayNames = parents.slice(0);

    var fixPathwayNames = function(pns){
      pathways = ["Pathways:"]
      overlaps = [" ", "Genes in Overlapping", "Pathways:"]
      abc = ["1: ", "2: ", "3: "]
      abcs = ["1 & 2", "1 & 3", "2 & 3", "1 & 2 & 3"]
      n = startingColors.length
      pns = pns.slice(0)
      check = 0
      for(i in pns){
        if(pns[i].includes(" & ")){
          check = 1
        }
      }
      if(!check){
        pns[pns.length-1] = "Novel Genes"
      }
      else{
        for(s in ordColors){
          color = ordColors[s]
          switch (color) {
            case color1:
              nS = abc[0] + pns[startingColors.indexOf(color)];
              pathways.push(nS)
              break;
            case color2:
              nS = abc[1] + pns[startingColors.indexOf(color)];
              pathways.push(nS)
              break;
            case color3:
              nS = abc[2] + pns[startingColors.indexOf(color)];
              pathways.push(nS)
              break;
            case novelColor:
              nS = "Novel Genes";
              pathways.push(nS)
              break;
            case color12:
              nS = abcs[0];
              overlaps.push(nS)
              break;
            case color13:
              nS = abcs[1];
              overlaps.push(nS)
              break;
            case color23:
              nS = abcs[2];
              overlaps.push(nS)
              break;
            case color123:
              nS = abcs[3];
              overlaps.push(nS)
              break;
          }
        }

        // bringing them together
        pnsSort = pathways.concat(overlaps)
        pns = pnsSort
      }
      return pns;
    }

    pathwayNames = fixPathwayNames(pathwayNames);

    legend.selectAll("text")
      .data(pathwayNames)
      .enter()
      .append("svg:text")
      .style("fill", function(d, i) {return legendColors[i];})
      .attr("id", function(d, i){return 'text' + i})
      .attr("transform", function(d) {return "translate(0,20)";})
      .attr("dy", function(d,i) {
        return 20*i
      })
      .text(function(d) { return d; });

    var path = svg.selectAll("path.link")
        .data(links)
      .enter().append("svg:path")
        .attr("class", function(d) { return "link source-" + d.source.key + " target-" + d.target.key; })
        .attr("d", function(d, i) { return line(splines[i]); })
        .style('stroke', function(d) {
          if(d.source.parent.name === 'Novel'){
            return d.target.color;
          }
          else{
            return d.source.color;
          }
        })
        .style("stroke-opacity", 0.05);

    var allNodes = svg.selectAll("g.node")
        .data(nodes.filter(function(n) { return !n.children; }))
      .enter().append("svg:g")
        .attr("class", "node2")
        .attr("id", function(d) { return "node2-" + d.key; })
        .attr("parentNode", function(d) { return d.parent.name;})
        .attr("transform", function(d) { return "rotate(" + (d.x - 90) + ")translate(" + d.y + ")"; });

    allNodes.append("svg:circle")
        .attr("r", Math.min(1/Math.log10(df.length)*5, 4))
        .style("fill", function(d) { return d.color; })
        .on("mouseover", mouseover)
        .on("mouseout", mouseout)
        .on("click", mouseclick);

    allNodes.append("svg:text")
        .style("fill", function(d){ return d.color;})
        .style("font-size", function(){
          return Math.min(1/Math.log10(df.length)*15, 14)
        })
        .attr("font-family", "HelveticaNeue-Light")
        .attr("dx", function(d) { return d.x < 180 ? 8 : -8; })
        .attr("dy", ".31em")
        .attr("text-anchor", function(d) { return d.x < 180? "start" : "end"; })
        .attr("transform", function(d) { return d.x < 180 ? null : "rotate(180)"; })
        .text(function(d) { return d.key; })
        .on("mouseover", mouseover)
        .on("mouseout", mouseout)
        .on("click", mouseclick);

    var still = false;

    svg.on("dblclick", mousedbl)
      .on("mousemove", mousemove)
      .on("mousedown", mousedown)
      .on("mouseup", mouseup);

    function mouse(e) {
      mouseX = e.pageX - rx - columnSpace;
      mouseY = e.pageY - ry - headSpace - niaidBannerY;
      return [mouseX, mouseY];
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

        rotate = rotate ? rotate+dm : dm;
        if (rotate > 360) rotate -= 360;
        else if (rotate < 0) rotate += 360;
        m0 = null;

        svg.style("-webkit-transform", null);

        svg
            .attr("transform", "translate(" + rx + "," + ry + ")rotate(" + rotate + ")")
          .selectAll("g.node text")
            .attr("dx", function(d) { return (d.x + rotate) % 360 < 180 ? 8 : -8; })
            .attr("text-anchor", function(d) { return (d.x + rotate) % 360 < 180 ? "start" : "end"; })
            .attr("transform", function(d) { return (d.x + rotate) % 360 < 180 ? null : "rotate(180)"; });
      }
    }

    function mouseover(d) {
      if(!still){
        if(!d.clicked){
          if(clickedData2.length === 0 || childrenArray2.indexOf(d.name) > -1){
            setvals(d, true, false)
          }
          if(clickedData2.length === 0){
            d3.selectAll(".vizText").remove()
            paintWindow(d, 1)
          }
          else if(childrenArray2.indexOf(d.name) > -1){
            d3.selectAll(".vizText").remove()
            paintWindow(d, 3);
          }
        }
        else{
            setvals(d, true);
        }
      }
    }

    function mouseout(d) {
      if(!still){
        var i = clickedData2.length>1 ? clickedData2.length-1 : 0;
        if(!d.clicked){
          setvals(d, false);
          if(clickedData2.length>=1){
            setvals(clickedData2[i], true, false)
          }
        }
        else if(clickedData2.length === 1){
          setvals(d, true, false)
          setvals(clickedData2[i], true, false)
        }
        if(clickedData2.length === 0){
          clearVizText()
        }
        else{
          d3.selectAll(".vizText").remove()
          paintWindow(d, 2)
        }
      }
    }

    function mouseclick(d){
      if(!still){
        d.clicked = !d.clicked && true;

        clickedData2.push(d)

        if(d.clicked){
          removelinks()

          setvals(d, true, false)

          d3.selectAll(".vizText").remove()
          paintWindow(d, 1);

          updateChildren(d)

          d3.selectAll(".node2 text").style("font-weight", "normal")

          childrenArray2 = childrenData2[d.name][0]

          d3.select("#node2-" + d.key + " text").style("font-weight", "bold")

          for(i in childrenArray2){
            childName = childrenArray2[i].split(".")[1]
            d3.select("#node2-" + childName + " text").style("font-weight", "bold")
          }
        }
        else if(clickedData2[clickedData2.length-1].name === d.name){
          setvals(d, true, false)
        }
        else{
          clickedData2.pop()

          updateChildren(d, true)

          setvals(d, false, false)
        }

        clickeRs = getClicker();

        Shiny.setInputValue("clickedData", clickeRs);
      }
    }

    function mousedbl(){
      childrenData2 = {}
      childrenArray2 = []
      still = false;
      rotate = null;

      svg.selectAll("path.link")
          .style("stroke-opacity", 0.05)

      if(clickedData2.length>0){
        for(var i in clickedData2){
          setvals(clickedData2[i], false, false)
          clickedData2[i].clicked = false
          svg.select("#node2-" + clickedData2[i].key)
              .style('font-weight', 'normal');
        }
        clickedData2 = []
      }

      svg.selectAll(".node2 text")
        .style("font-weight", "normal")
        .style("opacity", 1)

      svg.selectAll(".node2 circle").style('opacity', 1)

      svg.attr("transform",  "translate(" + rx + "," + ry + ")")
      clearVizText()

      Shiny.setInputValue("clickedData", '{}');
    }

    function clearVizText(){
      d3.selectAll(".vizText").remove()
      windowFields = ['Gene Name:', ' ',
                      'Interactions:',
                      'Screen Input:']
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
      childrenArray2 = []

      if(!remove){
        childrenData2 = {[d.name]: [d.imports]}
      }
      else{
        delete childrenData2[d.name]
      }
    }

    function getClicker(){
      nConn = clickedData2.length
      if(nConn===0){
        return '{}';
      }
      else if(nConn===1){
        clicker = clickedData2[0]
        clickeR = '{"Name1": ["' + clicker.key + '"],' +
                  '"Node1": ["' + clicker.name + '"],' +
                  '"Parent1": ["' + clicker.parent.name + '"],' +
                  '"Interactions1": ["' + clicker.imports.length + '"],' +
                  '"Confidence": ["' + clicker.Confidence + '"]}'
        return clickeR;
      }
      else{
        var clickeRs = '['
        for(conn in clickedData2){
            clicker = clickedData2[conn]
            if(conn < nConn-1){
              nextClicker = clickedData2[parseInt(conn)+1]
              clickeR = '{"Name1": ["' + clicker.key + '"],' +
                        '"Node1": ["' + clicker.name + '"],' +
                        '"Parent1": ["' + clicker.parent.name + '"],' +
                        '"Interactions1": ["' + clicker.imports.length + '"],' +
                        '"Confidence": ["' + clicker.Confidence + '"],' +
                        '"Name2": ["' + nextClicker.key + '"],' +
                        '"Node2": ["' + nextClicker.name + '"],' +
                        '"Parent2": ["' + nextClicker.parent.name + '"],' +
                        '"Interactions2": ["' + nextClicker.imports.length + '"],' +
                        '"Weight": ["' + clicker.weights[clicker.imports.indexOf(nextClicker.name)] + '"],' +
                        '"Source": ["' + clicker.datasource[clicker.imports.indexOf(nextClicker.name)] + '"]},'
            }
            else{
              clickeR = '{"Name1": ["' + clicker.key + '"],' +
                        '"Node1": ["' + clicker.name + '"],' +
                        '"Parent1": ["' + clicker.parent.name + '"],' +
                        '"Interactions1": ["' + clicker.imports.length + '"],' +
                        '"Confidence": ["' + clicker.Confidence + '"],' +
                        '"Name2": ["NA"], "Node2": ["NA"], "Parent2": ["NA"], "Interactions2": ["NA"], "Weight": ["NA"], "Source": ["NA"]}]'
            }
            clickeRs = clickeRs + clickeR
          }
        return clickeRs;
      }
    }

    function setvals(d, val, block=true, pop=8){
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

      svg.selectAll(".node2.source text")
          .filter(function(){
            return d3.select(this).attr("dx") == 8
          })
          .attr("dx", pop*val + pop)
          .style("font-weight", "bold")

      svg.selectAll(".node2.source text")
          .filter(function(){
            return d3.select(this).attr("dx") == -8
          })
          .attr("dx", -pop*val - pop)
          .style("font-weight", "bold")

      svg.selectAll(".node2")
          .filter(function(){
            return d3.select(this).select(" text").attr("dx") > 8 && d3.select(this).attr("class") === "node2";
          })
          .selectAll(" text")
          .attr("dx", 8)
          .style("font-weight", "normal")

      svg.selectAll(".node2")
          .filter(function(){
            return d3.select(this).select(" text").attr("dx") < -8 && d3.select(this).attr("class") === "node2";
          })
          .selectAll(" text")
          .attr("dx", -8)
          .style("font-weight", "normal")
    }

    function updateNodes(name, value) {
      return function(d) {
        if (value) this.parentNode.appendChild(this);

        lastclicked = clickedData2[clickedData2.length-1];
        //nodeBold = value || lastclicked===d[name] ? 'bold' : 'normal';

        svg.select("#node2-" + d[name].key)
          .classed(name, value);
      };
    }

    function removelinks(){
      for(var i = 0; i < clickedData2.length-1; i++){
        clickedData2[i].clicked = false

        svg.selectAll("path.link.target-" + clickedData2[i].key)
          .classed('target', false)
          .style("stroke-opacity", 0.05)
          .each(updateNodes("source", false));

        svg.selectAll("path.link.source-" + clickedData2[i].key)
          .classed('source', false)
          .style("stroke-opacity", 0.05)
          .each(updateNodes("target", false,));
      }
    }

    function revertOne(){
      if(still === false && clickedData2.length > 0){
        deleteGene = clickedData2[clickedData2.length-1]

        setvals(deleteGene, false, false)
        deleteGene.clicked = false

        d3.selectAll(".node2 text")
          .style("font-weight", "normal")

        d3.selectAll("path.link")
          .style("stroke-opacity", 0.05)

        clickedData2.pop()
        geneRevert = clickedData2[clickedData2.length-1]

        setvals(geneRevert, true, false)

        clearVizText()

        updateChildren(geneRevert)
        childrenArray2 = childrenData2[geneRevert.name][0]

        paintWindow(geneRevert, 1)
      }
    }

    function drawConnections(){
      still = true;
      nConn = clickedData2.length
      if(nConn === 0){
        still = false
        return;
      }
      else if(nConn === 1){
        svg.selectAll("path.link")
            .classed("source", false)
            .style("stroke-opacity", 0.01)

        svg.selectAll(".node2 text").style('opacity', 0.5)

        svg.selectAll(".node2 circle").style('opacity', 0.5)

        sourceName = clickedData2[0].key

        svg.selectAll("path.link.source-" + sourceName)
            .classed("source", true)
            .style("stroke-opacity", 1)

        svg.selectAll(".node2")
            .filter(function(){
              return d3.select(this).select(" text").attr("dx") != 8 &&
                      d3.select(this).select(" text").attr("dx") != -8 &&
                      d3.select(this).attr("class") !== "node";
            })
            .selectAll(" text")
            .style('opacity', 1)
            .style("font-weight", "bold")

        svg.selectAll(".node2")
            .filter(function(){
              return d3.select(this).select(" text").attr("dx") != 8 &&
                      d3.select(this).select(" text").attr("dx") != -8 &&
                      d3.select(this).attr("class") !== "node";
            })
            .selectAll(" circle")
            .style("opacity", 1)

      }
      else{
        svg.selectAll("path.link")
            .classed("source", false)
            .style("stroke-opacity", 0.01)

        svg.selectAll(".node2 text")
            .style("font-weight", "normal")
            .style('opacity', 0.5)

        svg.selectAll(".node2 circle").style('opacity', 0.5)

        svg.selectAll(".node2")
            .filter(function(){
              return d3.select(this).select(" text").attr("dx") < -8;
            })
            .selectAll(" text")
            .attr("dx", -8)

        svg.selectAll(".node2")
            .filter(function(){
              return d3.select(this).select(" text").attr("dx") > 8;
            })
            .selectAll(" text")
            .attr("dx", 8)

        for(i in clickedData2){
          i = parseInt(i)
          if(i != nConn-1){
            sourceName = clickedData2[i].key
            targetName = clickedData2[i+1].key
            svg.select("path.link.source-" + sourceName + ".target-" + targetName)
                .classed("source", true)
                .style("stroke-opacity", 1)
          }
          svg.select("#node2-" + clickedData2[i].key + " text")
              .attr("dx", function(d){
                if(d3.select("#node2-" + clickedData2[i].key + " text").attr("dx") == -8){
                  return -16;
                }
                else{
                  return 16;
                }
              })
              .style('font-weight', 'bold')
              .style('opacity', 1);

          svg.select("#node2-" + clickedData2[i].key + " circle")
              .style("opacity", 1)
        }
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
                        'Interactions: ' + d.datasource.length,
                        'Screen Input: ' + d.Confidence]
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
            .attr("font-family", "Helvetica")
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
        n = clickedData2[clickedData2.length-1]
        windowFields = ['Gene Name: ',
                        n.key,
                        'Interactions: ' + n.datasource.length,
                        'Screen Input: ' + n.Confidence]
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
            .attr("font-family", "Helvetica")
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
        n = clickedData2[clickedData2.length-1]
        cInd = childrenArray2.indexOf(d.name)
        windowFields = ['Linked Gene: ', d.key,
                        'Reference Gene: ', n.key,
                        'Interactions: ' + d.weights.length,
                        'Interaction Score: ' + n.weights[cInd],
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
            .attr("font-family", "Helvetica")
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
