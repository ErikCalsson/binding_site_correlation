# imports extern
import pandas as pd  # for data structuring
import dash  # visualisation of data via dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px

# imports intern
import src.data_calculation as data


# start dash
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

validated = 'Values'
validated += data.val_str  # result if p-value is bigger or smaller than chi²-value
validated += 'statistical significant different'
validated += '::'
validated += str(data.sci_out)  # TODO value for display only, remove later

# --------------------------------------------
# visualisation of output data

# dataframe for output
df = pd.DataFrame({
    "Overlap": ["Quotient", "First File", "Second File", "Quotient", "First File", "Second File"],
    "Coverage": [data.both_files_lazy, data.first_file_lazy, data.second_file_lazy,
                 data.both_files_log, data.first_file_log, data.second_file_log],
    "Size": ["Absolute", "Absolute", "Absolute",
             "Log 2", "Log 2", "Log 2"]
})

# figure
fig = px.bar(df, x="Overlap", y="Coverage", color="Size", barmode="group")

# app layout
app.layout = html.Div(children=[

    # title for the webpage
    html.H1(children="Overlap between both BED-files", style={'text-align': 'center'}),

    html.Div(id='output_container', children=[]),

    html.Br(),

    dcc.Graph(id='overlap_files', figure=fig),

    html.Br(),

    # html.H3(children=validated_own, style={'text-align': 'left'}),
    html.H3(children=validated, style={'text-align': 'left'})

])

# run app
app.run_server(debug=True)  # debug=true   means update browser by code change
# see in browser http://127.0.0.1:8050/ ONLY development server
