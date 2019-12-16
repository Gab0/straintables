#!/bin/python

import io
from flask import Response, Flask, request, render_template
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.pyplot import Figure
from straintables import Viewer, alignmentData, logo

import pathlib
import os
import argparse
import waitress

Description = ""
TemplateFolder = os.path.join(
    pathlib.Path(Viewer.__file__).resolve().parent,
    "WebComponents"
)
app = Flask("straintables Viewer", template_folder=TemplateFolder)


def GenerateLogoWeb():
    data = logo.genlogo()

    data = data.replace("<", "&lt;").replace(">", "&gt;")
    # .replace(" ", "&nbsp;")
    # return data.replace("\n", "<br>")
    return data


class ApplicationState():
    quad_allowed_regions = [0, 1]
    figsize = (15, 15)
    alnData = None
    Regions = None

    custom_allowed_regions = [0, 1]
    custom_reference_region = 0

    PlotOptions = {
        "fontsize": 12,
        "matrix_values": 0
    }

    def GetMatrixParameters(self):

        MatrixParameters = {
            "Normalize": False,
            "showNumbers": self.PlotOptions["matrix_values"]
        }
        return MatrixParameters


def parse_arguments():

    parser = argparse.ArgumentParser(description=Description)

    parser.add_argument('inputDir',
                        type=str,
                        nargs=1,
                        metavar="inputDirectory",
                        help='inputDirectory')

    parser.add_argument("-d",
                        metavar="inputDirectory",
                        dest="inputDirectory")

    parser.add_argument("--port",
                        type=int,
                        help="Local port to serve this http application.",
                        default=5000)
    options = parser.parse_args()

    if not options.inputDirectory:
        options.inputDirectory = options.inputDir[0]

    return options


@app.route("/alignment/<index>")
def ViewAlignment(index):
    pass


def BuildImage(fig, format="png") -> bytes:
    Formats = {
        "png": "print_png",
        "eps": "print_eps"
    }
    output = io.BytesIO()
    FigureCanvas(fig).__getattribute__(Formats[format])(output)

    return output.getvalue()


@app.route("/plot_custom/set")
def plot_custom_setup():
    AllowedRegionIndexes = [
        int(i)
        for i in request.args.getlist("allowed_regions")
    ]
    ReferenceIndex = [
        int(i)
        for i in request.args.getlist("reference_region")
    ]

    if ReferenceIndex:
        app.state.custom_allowed_regions = AllowedRegionIndexes
        app.state.custom_reference_region = ReferenceIndex[0]

    return plot_custom_view()


@app.route("/plot_custom/figure")
def plot_custom():
    fig = Figure(figsize=app.state.figsize)

    Viewer.plotViewport.plotRegionBatch(
        fig,
        app.state.alnData,
        app.state.custom_allowed_regions,
        reorganizeIndex=app.state.custom_reference_region,
        MatrixParameters=app.state.GetMatrixParameters())

    output = BuildImage(fig)
    return Response(output, mimetype='image/png')


@app.route("/plot_custom")
def plot_custom_view():
    return render_template("CustomPlotView.html",
                           logo=GenerateLogoWeb())


@app.route('/plot_quad/set')
def set_plot_quad():
    NewRegions = [request.args.get(a) for a in ["r1", "r2"]]

    print(NewRegions)
    for i, R in enumerate(NewRegions):
        if R is not None:
            app.state.quad_allowed_regions[i] = int(R)

    return render()


@app.route('/plot_quad/figure')
def plot_quad():

    fig = Figure(figsize=app.state.figsize)
    Viewer.plotViewport.MainDualRegionPlot(fig,
                                           app.state.alnData,
                                           app.state.quad_allowed_regions)
    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype='image/png')


@app.route("/debug")
def show_debug():
    message = str(app.state.custom_allowed_regions)
    message += str(app.state.PlotOptions["matrix_values"])
    return(Response(message))


@app.route("/")
def render():

    selected = app.state.quad_allowed_regions
    return render_template("MainView.html",
                           regions=app.state.Regions,
                           selected_regions=selected,
                           t=request.args.get("r1"),
                           logo=GenerateLogoWeb())


@app.route("/export")
def render_export():
    return render_template("CustomPlotBuild.html",
                           regions=app.state.Regions,
                           logo=GenerateLogoWeb(),
                           current_allowed_regions=app.state.custom_allowed_regions,
                           current_reference_region=app.state.custom_reference_region)


@app.route("/options")
def show_options_page():
    return render_template("PlotOptions.html",
                           options=app.state.PlotOptions,
                           logo=GenerateLogoWeb())

@app.route("/options/set")
def define_options():
    for opt in app.state.PlotOptions.keys():
        new_value = request.args.get(opt)
        print(new_value)
        if new_value is not None:
            app.state.PlotOptions[opt] = int(new_value)

    return render()


def main():
    app.state = ApplicationState()
    options = parse_arguments()
    app.state.alnData = alignmentData.AlignmentData(options.inputDirectory)
    app.state.Regions = app.state.alnData.MatchData["LocusName"]

    waitress.serve(app, port=options.port, url_scheme='http')
    #app.run(use_reloader=True, debug=True)


if __name__ == "__main__":
    main()
