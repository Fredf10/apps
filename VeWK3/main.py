'''
Created on Oct 5, 2016

@author: fredrik
'''



from classCreateSimpleApp import SetupApp

from bokeh.layouts import Row, column
from bokeh.io import curdoc


app = SetupApp()

for w in app.Widgetlist:
    w.on_change('value', app.update_data)


# Set up layouts and add to document
inputs = column(children=app.Widgetlist)
outputs = column(children=[app.plot_line])
curdoc().add_root(Row(children=[inputs, outputs], width=1200))