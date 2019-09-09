import os
import re

from flask import Blueprint, render_template, current_app, abort, make_response, redirect, url_for

import cea.inputlocator
import cea.plots
import cea.plots.categories

blueprint = Blueprint(
    'plots_blueprint',
    __name__,
    url_prefix='/plots',
    template_folder='templates',
    static_folder='static',
)

categories = {c.name: {'label': c.label, 'plots': [{'id': p.id(), 'name': p.name} for p in c.plots]}
              for c in cea.plots.categories.list_categories()}


@blueprint.route('/index')
def index():
    return redirect(url_for('plots_blueprint.route_react_dashboard'))


@blueprint.route('/dashboard')
def route_react_dashboard():
    debug = current_app.cea_config.get('general:debug')
    return render_template('react_dashboard.html', debug=debug, last_updated=dir_last_updated())


@blueprint.route('/dashboard/move_plot_up/<int:dashboard_index>/<int:plot_index>')
def route_move_plot_up(dashboard_index, plot_index):
    """Move a plot up in the dashboard"""
    if plot_index > 0:
        dashboards = cea.plots.read_dashboards(current_app.cea_config, current_app.plot_cache)
        dashboard = dashboards[dashboard_index]
        if plot_index < len(dashboard.plots):
            swap(dashboard.plots, plot_index - 1, plot_index)
            cea.plots.write_dashboards(current_app.cea_config, dashboards)
    return redirect(url_for('plots_blueprint.route_dashboard', dashboard_index=dashboard_index))


def swap(lst, i, j):
    """Swap positions of elements in a list as given by their indexes i and j"""
    lst[i], lst[j] = lst[j], lst[i]


@blueprint.route('/dashboard/move_plot_down/<int:dashboard_index>/<int:plot_index>')
def route_move_plot_down(dashboard_index, plot_index):
    """Move a plot down in the dashboard"""
    if plot_index >= 0:
        dashboards = cea.plots.read_dashboards(current_app.cea_config, current_app.plot_cache)
        dashboard = dashboards[dashboard_index]
        if plot_index < (len(dashboard.plots) - 1):
            swap(dashboard.plots, plot_index, plot_index + 1)
            cea.plots.write_dashboards(current_app.cea_config, dashboards)
    return redirect(url_for('plots_blueprint.route_dashboard', dashboard_index=dashboard_index))


@blueprint.route('/div/<int:dashboard_index>/<int:plot_index>')
def route_div(dashboard_index, plot_index):
    """Return the plot as a div to be used in an AJAX call"""
    try:
        plot = load_plot(dashboard_index, plot_index)
    except Exception as ex:
        return abort(500, ex)
    if not plot.missing_input_files():
        plot_div = plot.plot_div()
        # BUGFIX for (#2102 - Can't add the same plot twice in a dashboard)
        # update id of div to include dashboard_index and plot_index
        if plot_div.startswith("<div id="):
            div_id = re.match('<div id="([0-9a-f-]+)"', plot_div).group(1)
            plot_div = plot_div.replace(div_id, "{div_id}-{dashboard_index}-{plot_index}".format(
                div_id=div_id, dashboard_index=dashboard_index, plot_index=plot_index))
        return make_response(plot_div, 200)
    else:
        return render_template('missing_input_files.html',
                               missing_input_files=[lm(*args) for lm, args in plot.missing_input_files()],
                               script_suggestions=script_suggestions(
                                   lm.__name__ for lm, _ in plot.missing_input_files())), 404


@blueprint.route('/table/<int:dashboard_index>/<int:plot_index>')
def route_table(dashboard_index, plot_index):
    """Return the table for the plot as a div to be used in an AJAX call"""
    try:
        plot = load_plot(dashboard_index, plot_index)
    except Exception as ex:
        return abort(500, ex)
    if not plot.missing_input_files():
        return make_response(plot.table_div(), 200)


def script_suggestions(locator_names):
    """Return a list of CeaScript objects that produce the output for each locator name"""
    import cea.scripts
    schemas = cea.scripts.schemas()
    script_names = []
    for name in locator_names:
        script_names.extend(schemas[name]['created_by'])
    return [cea.scripts.by_name(n) for n in sorted(set(script_names))]


def load_plot(dashboard_index, plot_index):
    """Load a plot from the dashboard_yml"""
    cea_config = current_app.cea_config
    dashboards = cea.plots.read_dashboards(cea_config, current_app.plot_cache)
    dashboard_index = dashboards[dashboard_index]
    plot = dashboard_index.plots[plot_index]
    return plot


def dir_last_updated():
    return str(max(os.path.getmtime(os.path.join(root_path, f))
                   for root_path, dirs, files in os.walk(os.path.join(os.path.dirname(__file__), 'static'))
                   for f in files))
