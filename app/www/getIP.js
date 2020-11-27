$(document).ready(function(){
  $.get("https://ipinfo.io", function(response) {
    Shiny.onInputChange("getIP", response);
  }, "json");
});
