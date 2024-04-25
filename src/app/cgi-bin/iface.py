#!/var/www/wrapper.sh
#url: http://localhost/app/cgi-bin/iface.py?method=runSimulation&start_date=2017-7-1&end_date=2018-1-1&location=cordoba&scenario=rc
import cgi
import spaa

#main
#parse parameters
params = {}
fieldStorage=cgi.FieldStorage()
for key in fieldStorage.keys():
   params[ key ] = fieldStorage[ key ].value
#make a callable method
method=getattr(spaa,params['method'])

#call the method
print('Content-Type: application/json\n\n')
print(method(params))
