"use strict";
var app={};
$(document).ready(function() {

  app.init= function(){
    $("#runSimulation").click(app.runSimulation);
    //Configure the datepicker
    $(".date").datepicker({
      format: 'yyyy-mm-dd',
      endDate: '+30d',
      autoclose: true,
      keyboardNavigation:false,
    });
    //Disable submit forms via enter key
    $('form').bind('keydown', function(e) {
      if (e.keyCode == 13) {
          e.preventDefault();
      }
    });

    $('#start_date').datepicker('update','2015-8-1');
    var end_date=new Date();
    end_date.setDate(end_date.getDate() + 29);//Doesn't take the last day.conflict with validation endDate: '+30d'.
    $("#end_date").datepicker("update",end_date);
  };

  app.runSimulation=function(e){
    e.preventDefault();
    $('.alert').hide();
    $('#population').html('Loading...');
    $('#population_II').html('Loading...');
    $('#weather').html('Loading...');
    var parameters=new Object();
    parameters.start_date=$('#start_date').val();
    parameters.end_date=$('#end_date').val();
    parameters.location=$('#location').val();
    parameters.scenario=$('#scenario').val();
    $.ajax({
      type: "GET",
      contentType: "application/json",
      url: "cgi-bin/iface.py?method=runSimulation&"+$.param( parameters ),
      success: function(data) {
        app.draw_chart('population',data.population);
        app.draw_chart('weather',data.weather);
        app.draw_chart('population_II',data.population_II);
        $('#download').attr('href','cgi-bin/iface.py?method=download&filename='+data.filename);
      },
      error:function(data){
        $('.alert #message').html('Error inesperado.'+data.responseText.substring(0,500));
        $('.alert').show();
        console.log(data.responseText);
      }
    });
  }


  app.draw_chart=function(plot_container_id,series_xy) {
    console.log('data.population_II',series_xy);
    var traces=[];
    $(series_xy).each(function(index,serie_xy){
      var trace={
        x:[],
        y:[],
        name: serie_xy.name,
        type: serie_xy.type
      }
      $(serie_xy.data).each(function( index, value ) {
        trace.x.push(new Date(value[0]));
        trace.y.push(value[1]);
      });
      traces.push(trace);
    });
    traces.push({x: [new Date()],y: [0],mode: 'markers',name: 'Fecha Actual',type: 'scatter'});//today marker

    $('#'+plot_container_id).html('');
    Plotly.newPlot(plot_container_id, traces);
  }

  app.init();
});
