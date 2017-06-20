'''
Created on Oct 3, 2016

@author: fredrik
'''

from classSetupDirApp import SetupApp

from bokeh.io import curdoc
from bokeh.layouts import Row, column

app = SetupApp()

for w in app.Widgetlist:
    w.on_change('value', app.update_data)



# Set up layouts and add to document
inputs = column(children=app.Widgetlist)
outputs = column(children=[app.plotP, app.plotQ])
outputs2 = column(children=[app.plot_Img])
curdoc().add_root(Row(children=[inputs, outputs, outputs2], width=2200))