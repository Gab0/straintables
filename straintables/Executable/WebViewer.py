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


def GenerateBlankFigure(PlotOptions):
    figsize = PlotOptions["figsize"]
    return Figure(figsize=(figsize, figsize), dpi=PlotOptions["dpi"])


def GenerateLogoWeb():
    data = logo.genlogo()

    data = data.replace("<", "&lt;").replace(">", "&gt;")
    # .replace(" ", "&nbsp;")
    # return data.replace("\n", "<br>")
    return data


class ApplicationState():
    # Global Options
    alnData = None
    Regions = None

    # Global Configurable Options
    PlotOptions = {
        "fontsize": 12,
        "dpi": 100,
        "matrix_values": 0,
        "figsize": 15,
        "format": "png",
        "YlabelsOnRight": False
    }

    # Quad mode options
    quad_allowed_regions = [0, 1]

    # Custom mode options
    custom_allowed_regions = [0, 1]
    custom_reference_region = 0

    def GetMatrixParameters(self):

        MatrixParameters = {
            "Normalize": False,
            "showNumbers": self.PlotOptions["matrix_values"],
        }

        MatrixParameters.update(self.PlotOptions)
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
                        dest="inputDirectory",
                        help="Input directory with analysis data.")

    parser.add_argument("--port",
                        type=int,
                        help="Local port to serve this http application.",
                        default=5000)

    parser.add_argument("--debug",
                        dest="Debug",
                        action="store_true")

    options = parser.parse_args()

    if not options.inputDirectory:
        options.inputDirectory = options.inputDir[0]

    return options


@app.route("/alignment/<index>")
def ViewAlignment(index):
    pass


def BuildImage(fig, imgFormat="png") -> Response:
    Formats = {
        "png": "print_png",
        "eps": "print_eps"
    }
    output = io.BytesIO()

    for i in range(3):
        fig.tight_layout()
    FigureCanvas(fig).__getattribute__(Formats[imgFormat])(output)

    return Response(output.getvalue(), mimetype="image/%s" % imgFormat)


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
    fig = GenerateBlankFigure(app.state.PlotOptions)
    app.state.CurrentPlots = Viewer.plotViewport.plotRegionBatch(
        fig,
        app.state.alnData,
        app.state.custom_allowed_regions,
        reorganizeIndex=app.state.custom_reference_region,
        MatrixParameters=app.state.GetMatrixParameters())

    return BuildImage(fig, app.state.PlotOptions["format"])


@app.route("/plot_custom")
def plot_custom_view():
    return render_template("CustomPlotView.html",
                           logo=GenerateLogoWeb())


# -- Plot quad methods;
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
    fig = GenerateBlankFigure(app.state.PlotOptions)

    app.state.CurrentPlots = Viewer.plotViewport.MainDualRegionPlot(
        fig,
        app.state.alnData,
        app.state.quad_allowed_regions,
        MatrixParameters=app.state.GetMatrixParameters()
    )
    app.state.CurrentFigure = fig
    return BuildImage(fig, app.state.PlotOptions["format"])


@app.route('/export_eps')
def export_eps():
    Path = request.args.get("export_directory")
    filepaths = []

    # -- Export main figure;
    if len(app.state.CurrentFigure) > 1:
        filename = "figure_full.eps"
        output_path = os.path.join(Path, filename)
        app.state.CurrentFigure.savefig(output_path)
        filepaths.append(filename)

    # -- Export individual matrices;
    for p, plotableMatrix in enumerate(app.state.CurrentPlots):
        fig = GenerateBlankFigure(app.state.PlotOptions)
        Viewer.plotViewport.createMatrixSubplot(fig, 111, plotableMatrix)

        filename = "figure%i_%s.eps" % (p + 1, plotableMatrix.name)
        output_path = os.path.join(Path, filename)
        print("Writing %s." % output_path)

        filepaths.append(output_path)
        fig.savefig(output_path)

    m = "Sucessfully created the following files:"
    Message = "\n<br>".join([m] + filepaths)

    return(Response(Message))


@app.route("/debug")
def show_debug():
    def show_dict(data):
        m = ""
        for k in data.keys():
            m += "%s: %s<br>\n" % (k, data[k])
        return m

    message = "<html>\n\n"
    message += str(app.state.custom_allowed_regions)
    message += str(app.state.PlotOptions["matrix_values"])

    message += show_dict(app.state.__dict__)
    matrix_options = app.state.GetMatrixParameters()
    message += show_dict(matrix_options)

    message += "</html>"
    return(Response(message))


@app.route("/")
def render():

    selected = app.state.quad_allowed_regions

    region_names = app.state.alnData.getRegionNamesFromIndex(app.state.quad_allowed_regions)
    currentPWMData = app.state.alnData.findPWMDataRow(*region_names)
    MatchData = app.state.alnData.getMatchDataFromNames(region_names)
    information = Viewer.plotViewport.RegionData(currentPWMData,
                                                 MatchData, *region_names)
    information = "<br>".join(information)
    return render_template("MainView.html",
                           information=information,
                           regions=app.state.Regions,
                           selected_regions=selected,
                           t=request.args.get("r1"),
                           logo=GenerateLogoWeb())


@app.route("/export")
def render_export():
    return render_template(
        "CustomPlotBuild.html",
        regions=app.state.Regions,
        logo=GenerateLogoWeb(),
        current_allowed_regions=app.state.custom_allowed_regions,
        current_reference_region=app.state.custom_reference_region
    )


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
            try:
                V = int(new_value)
            except ValueError:
                V = new_value
            app.state.PlotOptions[opt] = V

    return render()


def Execute(options):
    app.state = ApplicationState()
    app.state.alnData = alignmentData.AlignmentData(options.inputDirectory)
    app.state.Regions = app.state.alnData.MatchData["LocusName"]

    print(logo.logo)
    print()
    print("loading straintables matrix viewer server...")
    print("\t point your web browser to the address below.")
    print()

    if options.Debug:
        app.run(use_reloader=True, debug=True)
    else:
        waitress.serve(app, port=options.port, url_scheme='http')


def main():
    options = parse_arguments()
    Execute(options)


if __name__ == "__main__":
    main()
