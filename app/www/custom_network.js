Shiny.addCustomMessageHandler("jsondata",
  function(message){
    var json_data = message;
    console.log(json_data)
  });
