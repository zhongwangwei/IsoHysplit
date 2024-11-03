from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import os
import matplotlib.ticker as tk
import matplotlib.colors as clr

def map_labeller(basemap, which, labels, labelstyle, labelzorder):
    """
    Write labels on ``basemap``.

    Use ``labelfile_generator()`` and ``labelfile_reader()`` to get
    ``labels`` and ``labelstyle`` in the appropriate order and format.

    Parameters
    ----------
    basemap : ``Basemap`` instance
        Any ``Basemap`` instance.  For easy map creation, see ``MapDesign``
        class.
    which : list of strings
        Can include ['sea'|'country'|'ocean'|'place'|'city'].
        The labels to apply
    labels : list of lists of strings or tuples of floats
        List of lists of labels and coordinates
    labelzorder : int
        Zorder of map labels.
    labelstyle : list of dictionaries
        List of dictonaries that contain label formatting parameters

    """
    # Initialize lists of labels and coordinates
    (sea_labels, sea_coords,
     cities_labels, cities_coords,
     country_labels, country_coords,
     ocean_labels, ocean_coords,
     place_labels, place_coords) = labels

    (sea_style, cities_style, country_style, ocean_style,
     place_style) = labelstyle

    # Initialize formatting dictionaries
    sea_dict = {'place': sea_labels,
                'coord': sea_coords,
                'ha': 'center',
                'va': 'center',
                'fs': sea_style['fontsize'],
                'wt': sea_style['weight'],
                'fst': sea_style['fontstyle'],
                'zrd': labelzorder,
                'off': (0, 0),
                'add': ''}

    ocean_dict = {'place': ocean_labels,
                  'coord': ocean_coords,
                  'ha': 'center',
                  'va': 'center',
                  'fs': ocean_style['fontsize'],
                  'wt': ocean_style['weight'],
                  'fst': ocean_style['fontstyle'],
                  'zrd': labelzorder,
                  'off': (0, 0),
                  'add': ''}

    country_dict = {'place': country_labels,
                    'coord': country_coords,
                    'ha': 'center',
                    'va': 'center',
                    'fs': country_style['fontsize'],
                    'wt': country_style['weight'],
                    'fst': country_style['fontstyle'],
                    'zrd': labelzorder,
                    'off': (0, 0),
                    'add': ''}

    place_dict = {'place': place_labels,
                  'coord': place_coords,
                  'ha': 'right',
                  'va': 'center',
                  'fs': place_style['fontsize'],
                  'wt': place_style['weight'],
                  'fst': place_style['fontstyle'],
                  'zrd': labelzorder,
                  'off': (0.0, 1.0),
                  'add': r'$\bigstar$'}

    city_dict = {'place': cities_labels,
                 'coord': cities_coords,
                 'ha': 'right',
                 'va': 'center',
                 'fs': cities_style['fontsize'],
                 'wt': cities_style['weight'],
                 'fst': cities_style['fontstyle'],
                 'zrd': labelzorder,
                 'off': (0.0, 1.0),
                 'add': r'$\bullet$'}

    label_dict = {'sea': sea_dict,
                  'city': city_dict,
                  'place': place_dict,
                  'ocean': ocean_dict,
                  'country': country_dict}

    maplabels = [label_dict[label] for label in which]

    # Plot labels
    for j in maplabels:

        for i in range(0, len(j['place'])):

            x, y = basemap(j['coord'][i][1] + j['off'][1],
                           j['coord'][i][0] + j['off'][0])

            basemap.ax.text(x, y, j['place'][i] + j['add'],
                            horizontalalignment=j['ha'],
                            verticalalignment=j['va'],
                            fontsize=j['fs'], weight=j['wt'],
                            fontstyle=j['fst'], zorder=j['zrd'])


def labelfile_reader(labelfile):
    """
    Open and read a text file with label information.

    Parameters
    ----------
    labelfile : string
        Full or relative path to text file containing label information

    Returns
    -------
    labels : list of lists
        List of lists of labels and coordinates
    labelstyle : list of dictionaries
        List of dictonaries that contain label formatting parameters

    """
    labels = []
    labelstyle = []
    labeltypes = ['SEA', 'CITY', 'COUNTRY', 'PLACE', 'OCEAN']

    # Look for file
    if not os.path.exists(labelfile):
        raise OSError(labelfile, ' does not exist.\n',
                      'Generate a template at this location',
                      'using labelfile_generator(', labelfile, ')')

    # Open file
    with open(labelfile, 'r') as lfile:

        # Read first line, a label header
        line = lfile.readline().strip()

        # Cycle through types of labels, breaking at end of file
        while True:
            if line == '':
                break
            labellist = []
            coordlist = []

            # Run through labels within type
            while True:
                line = lfile.readline().strip()

                # If line is a new type
                if line in labeltypes:
                    # Get the next line, the style parameters
                    line = lfile.readline().strip()

                    if line == '':
                        break

                    # Break into substrings, one parameter per string
                    line = line.split('  ', 2)

                    # Split again over the '=', arrange into dictionary
                    for i in range(0, len(line)):
                        line[i] = line[i].strip().split('=')

                    styledict = {k: v for k, v in line}
                    labelstyle.append(styledict)

                    break

                elif line == '':
                    break

                else:
                    latitude = float(line[:6])
                    longitude = float(line[6:14])

                    coord = (latitude, longitude)

                    place = line[19:]
                    place = place.replace('\\n', '\n')

                    labellist.append(place)
                    coordlist.append(coord)

            if len(labellist) > 0:
                labels.append(labellist)
                labels.append(coordlist)

    return labels, labelstyle


def labelfile_generator(labelfile, example='east'):
    """
    Generate a label file template.

    Parameters
    ----------
    labelfile : string
        Full or relative path to template location
    example : string
        Default 'east'.  Also accepts 'west'.  Sample
        label file for each hemisphere.  Each has slightly different
        formatting.

    """
    # Open new file
    lbdict = {'east': ['File header\n',
                       'SEA\n',
                       '  fontstyle=italic   weight=normal   fontsize=16\n',
                       '  16.00    88.50     Bay of\\nBengal\n',
                       ' -15.05   115.00     South\\nChina\\nSea\n',
                       '  27.00   125.00     East\\nChina\\nSea\n',
                       'CITY\n',
                       '  fontstyle=normal   weight=normal   fontsize=15\n',
                       '  39.91   116.39     Beijing\n',
                       '  32.05   118.77     Nanjing\n',
                       '  25.27   110.28     Guilin\n',
                       'COUNTRY\n',
                       '  fontstyle=normal   weight=bold     fontsize=20\n',
                       '  35.00   100.00     CHINA\n',
                       'OCEAN\n',
                       '  fontstyle=italic   weight=bold     fontsize=20\n',
                       '  27.00   150.00     Pacific\\n\\nOcean\n',
                       '  -5.00    70.00     Indian Ocean\n',
                       'PLACE\n',
                       '  fontstyle=normal   weight=normal   fontsize=15\n',
                       '  32.29   119.05     Hulu Cave\n'],
              'west': ['File header\n',
                       'SEA\n',
                       '  fontstyle=italic   weight=normal   fontsize=13\n',
                       '  26.00   -90.00     Gulf of\\nMexico\n',
                       'CITY\n',
                       '  fontstyle=normal   weight=normal   fontsize=15\n',
                       '  44.98   -93.26     Minneapolis\n',
                       'COUNTRY\n',
                       '  fontstyle=normal   weight=normal   fontsize=15\n',
                       '  40.00   -110.00    United\\nStates\n',
                       'OCEAN\n',
                       '  fontstyle=italic   weight=normal   fontsize=18\n',
                       '  25.00  -150.00     Pacific\\n\\nOcean\n',
                       'PLACE\n',
                       '  fontstyle=normal   weight=normal   fontsize=15\n',
                       '  38.98  -114.30     Great Basin\\nNational Park\n']}

    if example is not 'east' and example is not 'west':
        example = 'east'

    labels = lbdict['example']

    with open(labelfile, 'w') as labelfile:

        labelfile.writelines(labels)
        labelfile.flush()


def traj_scatter(data, lons, lats, hymap, zorder=19, colormap='default',
                 edgecolor='none', size=25, sizedata=None, cnormalize=None,
                 snormalize=None, vmin=None, vmax=None, levels=11,
                 suppress_printmsg=False, **kwargs):
    """
    Scatter-plot of ``Trajectory``, ``TrajectoryGroup``, or ``Cluster`` data.
    
    Parameters remain the same, except:
    - hymap is now expected to be a Cartopy GeoAxes instance instead of Basemap
    - removed latlon parameter as it's no longer needed with Cartopy
    """
    norm = None
    msg = ('Use `cbar.ax.set_yticklabels()` ' +
           'or cbar.ax.set_xticklabels()` to change tick labels')

    transform_dict = {'sqrt': np.sqrt,
                      'log': np.log10,
                      'ln': np.log}

    if colormap is 'default':
        try:
            colormap = plt.cm.viridis
        except AttributeError:
            colormap = plt.cm.jet


    if cnormalize is 'boundary':
        if vmin is None:
            vmin = data.min()
        if vmax is None:
            vmax = data.max()
        bounds = np.linspace(vmin, vmax, levels)
        norm = clr.BoundaryNorm(bounds, colormap.N)
    elif cnormalize is 'log':
        norm = clr.LogNorm(vmin=vmin, vmax=vmax)
    elif cnormalize is 'ln':
        data = np.log(data)
        if not suppress_printmsg:
            print(msg, '\nnatural log normalization')
    elif cnormalize is 'sqrt':
        data = np.sqrt(data)
        if not suppress_printmsg:
            print(msg, '\nsqrt normalization')
    elif cnormalize is not None:
        try:
            norm = clr.PowerNorm(cnormalize, vmin=vmin, vmax=vmax)
        except:
            pass

    if sizedata is not None:
        if snormalize is not None:
            sizedata = transform_dict[snormalize](sizedata)
        size = sizedata * size

    # Use Cartopy's transform for plotting coordinates
    cm = hymap.scatter(lons, lats, c=data, s=size, cmap=colormap,
                      vmin=vmin, vmax=vmax, zorder=zorder,
                      edgecolor=edgecolor, norm=norm,
                      transform=ccrs.PlateCarree(), **kwargs)

    return cm


def meteo_contouring(hymap, data, longitudes, latitudes, contourf=True,
                     vmin=None, vmax=None, steps=50, levels=None, colors=None,
                     colormap=plt.cm.nipy_spectral, zorder=13, **kwargs):
    """
    Create contour or filled contour maps of ``data``.

    Parameters
    ----------
    hymap : ``Basemap`` instance
        Any ``Basemap`` instance.  For easy map creation, see ``MapDesign``
        class
    data : (M, N) ndarray of floats
        The information to contour
    longitudes : (M) ndarray of floats
        X-coordinates of ``data`` in decimal degrees
    latitudes : (N) ndarray of floats
        Y-coordinates of ``data`` in decimal degrees
    contourf : Boolean
        Default ``True``.  Create filled contour (``True``) or contour
        (``False``) plot.
    vmin : int or float
        Default ``None``.  The minimum value for contouring.
        If ``None``, ``vmin`` is the ``data`` minimum.
    vmax : int or float
        Default ``None``.  The maximum value for contouring.
        If ``None``, ``vmax`` is the ``data`` maximum.
    steps : int
        Default 50.  The number of steps between ``vmin`` and ``vmax``.
    levels : list of ints or floats
        Default ``None``.  The contouring levels, overriding level creation
        with ``vmin``, ``vmax``, ``steps``
    colors : list of strings or tuples
        Default ``None``.  The colors to use for contouring.
    colormap : ``matplotlib`` colormap
        Default ``plt.cm.Blues``.  Any ``matplotlib`` colormap.
    zorder : int
        Default 13.  Zorder of ``data`` on ``hymap``.
    **kwargs
        Passed to ``Basemap.contour()`` then ``Axes.contour()``
        (or ``Axes.contourf()``)

    Returns
    -------
    cm : ``matplotlib.contour.QuadContourSet`` instance
        Mappable for use in creating colorbars.  Colorbars may be created
        in ``PySPLIT`` using ``make_cbar()`` or ``make_cax_cbar()``

    """
    if vmin is None:
        vmin = data.min()
    if vmax is None:
        vmax = data.max()

    if levels is None:
        levels = np.linspace(vmin, vmax, steps)

    if colors is not None:
        colormap = None

    if longitudes.ndim == 1 and data.ndim == 2:
        longitudes, latitudes = np.meshgrid(longitudes, latitudes)

    if contourf:
        cm = hymap.contourf(longitudes, latitudes, data, zorder=zorder,
                            cmap=colormap, levels=levels, latlon=True,
                            colors=colors, **kwargs)
    else:
        cm = hymap.contour(longitudes, latitudes, data, zorder=zorder,
                           levels=levels, colors=colors, latlon=True,
                           cmap=colormap, **kwargs)

    return cm

def adjust_contourparams(cm, contours, colors=[None],
                         othercontours_visible=True, **kwargs):
    """
    Adjust contour parameters.

    Shortcut for recoloring particular contours and/or rendering other
    contours invisible.  Can also pass other kwargs to chosen contours.

    Parameters
    ----------
    cm : ``matplotlib.contour.QuadContourSet`` instance
        The contour set to adjust
    contours : list of ints or floats
        The levels to adjust
    colors : list of strings, tuples
        Default [``None``].  The colors of ``contours``
    othercontours_visible : Boolean
        Default ``True``.  If ``False``, then levels not in ``contours`` will
        be set invisible.
    **kwargs
        Collection of keywords for ``contours``.

    """
    if len(contours) != len(colors):
        colors = [colors[0]] * len(contours)

    for level, coll in zip(cm.levels, cm.collections):
        if level in contours:
            ind = contours.index(level)
            color = colors[ind]
            plt.setp(coll, **kwargs)
            if color is not None:
                coll.set_color(color)
        else:
            if not othercontours_visible:
                coll.set_alpha(0)

def make_cbar(data, ax, orientation='horizontal', cbar_size=(20, 1.0),
              reverse_cbar=False, **kwargs):
    """
    Make a colorbar on the same axis as ``ax``.

    Parameters
    ----------
    data : ``matplotlib PathCollection``
        The mappable
    ax : ``Axes`` instance
        The axis on which ``data`` is plotted
    orientation : string
        Default 'horizontal'.  ['horizontal'|'vertical'].  Colorbar orientation
    cbar_size : tuple of int, float
        Default (20, 1.0).  Colobar (aspect, shrink).  The H/W ratio of the
        colorbar, the fractional size of the colorbar.
    reverse_cbar : Boolean
        Default ``False``. If ``True``, colorbar is flipped over short axis.
        Value-color mapping is unaffected.
    **kwargs
        Passed to ``edit_cbar()``

    Returns
    -------
    cbar : ``matplotlib`` ``ColorBar`` instance
        The new colorbar

    """
    # Initialize colorbar
    cbar = plt.colorbar(data, ax=ax, orientation=orientation,
                        aspect=cbar_size[0], shrink=cbar_size[1])

    # Reverse colorbar
    if reverse_cbar:
        if orientation is 'horizontal':
            cbar.ax.invert_xaxis()
        else:
            cbar.ax.invert_yaxis()

    # Make pretty
    edit_cbar(cbar, **kwargs)

    return cbar

def make_cax_cbar(fig, rect, data, orientation='horizontal',
                  reverse_cbar=False, extend='neither', **kwargs):
    """
    Make a colorbar on a new axis.

    Parameters
    ----------
    fig : ``figure`` instance
    rect : list of floats
        The colorbar position and size.  [Distance from left, distance from
        bottom, size in x dimension, size in y dimension]
    data : ``matplotlib PathCollection``
        Mappable
    orientation : string
        Default 'horizontal'.  ['horizontal'|'vertical'].  The orientation of
        the colormapping within the colorbar.
    cbar_size : tuple of int, float
        Default (20, 1.0).  Colobar (aspect, shrink).  The H/W ratio of the
        colorbar, the fractional size of the colorbar.
    reverse_cbar : Boolean
        Default ``False``. If ``True``, colorbar is flipped over short axis.
        Value-color mapping is unaffected.
    extend : string
        Default 'neither'.  ['both'|'neither'|'under'|'over'].
        Extend colorbar with pointed ends.
    **kwargs
        Passed to ``edit_cbar()``

    Returns
    -------
    cax : ``matplotlib Axes`` instance
        The axis of the new colorbar.  Remove using ``fig.delaxes(cax)``
    cbar : ``matplotlib ColorBar`` instance
        The new colorbar

    """
    # Initialize cax and colorbar on cax
    cax = fig.add_axes(rect)
    cbar = fig.colorbar(data, cax=cax, orientation=orientation,
                        extend=extend)
    # Reverse colorbar
    if reverse_cbar:
        if orientation is 'horizontal':
            cbar.ax.invert_xaxis()
        else:
            cbar.ax.invert_yaxis()

    # Make pretty
    edit_cbar(cbar, **kwargs)

    return cax, cbar

def edit_cbar(cbar, divisions=5, cbar_label=None, tick_fs=16, label_fs=18,
              labelpad=24, rotation=0, tick_dir='out', tick_dim=(4, 2)):
    """
    Make the colorbar pretty.

    Adjust fontsizes, add label, get a reasonable number of nice ticks, etc.

    Parameters
    ----------
    cbar : ``matplotlib colorbar`` instance
        The colorbar created in ``make_cbar()`` or ``make_cax_cbar()``.
    divisions : int
        Default 5.  The number of nice ticks on the colorbar.  May be ``None``.
    cbar_label : string
        Default ``None``.  Colorbar label.
    tick_fs : int
        Default 16.  Font size of ticks
    label_fs : int
        Default 18.  Font size of ``cbar_label``
    labelpad : int
        Default 24.  Spacing between tick labels and ``cbar`` label
    rotation : int
        Default 0.  Label rotation in degrees.
    tick_dir : string
        Default 'out'.  ['out'|'in'|'inout']
        Direction that ticks are pointing relative to colorbar
    tick_dim : tuple of floats
        Default (4, 2).  The (length, width) of ``cbar`` ticks

    """
    # Adjust ticks and tick labels
    if divisions is not None:
        cbar.locator = tk.MaxNLocator(divisions, integer=False)

    cbar.ax.tick_params(labelsize=tick_fs, direction=tick_dir,
                        length=tick_dim[0], width=tick_dim[1])
    cbar.update_ticks()

    # Label colorbar
    if cbar_label is not None:
        cbar.set_label(cbar_label, labelpad=labelpad, fontsize=label_fs,
                       rotation=rotation)

    # Cbar will have lines through it if mappable's alpha < 1
    cbar.set_alpha(1)
    try:
        cbar.draw_all()
    except AttributeError:
        try:
            cbar._draw_all()  # For newer matplotlib versions
        except AttributeError:
            pass  # If neither method exists, skip drawing

def random_colors(number_ofcolors):
    """
    Generate random RGB tuples.

    Parameters
    ----------
    number_ofcolors : int
        Number of tuples to generate

    Returns
    -------
    colors : list of tuples of floats
        List of ``len(number_ofcolors)``, the requested random colors
    """
    color_tmp = np.random.rand(number_ofcolors, 3)
    color_tmp = np.vsplit(color_tmp, number_ofcolors)
    colors = []
    for c in color_tmp:
        colors.append(c[0])

    return colors

class MapDesign(object):
    """Class for holding map design elements."""

    def __init__(self, mapcorners, standard_pm, projection='PlateCarree',
                 mapcolor='light', maplabels=None, area_threshold=10000,
                 resolution='110m', zborder=14, zlandfill=12,
                 zmapbound=10, zlatlon=11, lat_labels=['left'],
                 lon_labels=['top'], latlon_labelspacing=(10, 20),
                 latlon_fs=20, latlon_spacing=(10, 20), drawstates=False,
                 drawoutlines=True, draw_latlons=True, land_alpha=0.85):
        """
        Initialize ``MapDesign`` instance.

        Is your map blank?  Try changing zmapbound.

        Parameters
        ----------
        mapcorners : list of floats
            Used to construct the map view for 'conic' and 'cyl' projections.
            Lower left longitude, latitude; upper right longitude, latitude.
        standard_pm : list of floats
            For cylindrical and conic projections, the list creates standard
                parallels and meridians
                (``lon_0``, ``lat_0``, ``lat_1``, ``lat_2``).
            For orthographic projection, ``lon_0`` and ``lat_0`` only
                are required.  Sets the view from above the earth.
            For polar projections, ``lon_0`` indicates the longitude that
                will be oriented N-S. ``lat_0`` is replaced by
                the ``boundinglat``, the lowest latitude that should appear
                on the map.  ``lat_1`` and ``lat_2`` not required
        projection : string
            Indicates map projection.  Default 'PlateCarree'.
                'cyl' : Equidistant cylindrical
                'cea' : Equal Area cylindrical
                'lcc' : Lambert Conformal Conic
                'aea' : Albers Equal Area Conic
                'ortho' : Orthographic (globe)
                'npstere' : North polar steroegraphic (conformal)
                'spstere' : South polar steroegraphic (conformal)
                'nplaea' : North polar azimuthal (equal area)
                'splaea' : South polar azimuthal (equal area)
        mapcolor : string
            Default 'light'. The map grayscheme.
            ['light'|'medium'|'dark'|None]
            Not available for 'ortho' projections
        maplabels : tuple of strings
            Default ``None``.
            (Label group, label file full/relative path, optional: zorder).
            Label group is a list of any or all of:
            ['sea', 'city', 'country','ocean','place']
            Label zorder defaults to 15.
        area_threshold : int
            Default 10000.  The minimum surface area a feature must have to
            be drawn on the map.
        resolution : char
            Default 'c'.  ['c'|'l'|'i'|'h'|'f'].
            Crude, low, intermediate, high, full. The relative resolution of
            map boundaries.  Drops off by about 80 percent between datasets.
        zborder : int
            Default 14. The zorder of country and coastal outlines.
        zlandfill : int
            Default 12.  The zorder of the continent fill color.
        zmapbound : int
            Default 16.  The zorder of the map boundary in older versions of
            basemap (background correctly put at bottom of stack),
            and zorder of the ocean/background in newer versions (boundary
            correctly put at top of stack).
            Try zmapbound=10 if 16 yields a blank map.
        zlatlon : int
            Default 11. The zorder of the lines of latitutde and longitude.
        lat_labels : list of strings
            Default ['right'].
            The sides of the map that should have the latitudes labelled.
        lon_labels : list of strings
            Default ['top'].
            The sides of the map that should have the longitudes labelled.
        latlon_labelspacing : int or float
            Default (10, 20).  Degrees between (latitude, longitude) labels
        latlon_fs : int or float
            Default 20.  Font size of latitude, longitude labels.
        latlonspacing : int or float
            Default (10, 20).  Degrees between plotted lines of latitude.
        drawstates : Boolean
            Default False.  Draw state outlines on ``Basemap``.
        drawoutlines : Boolean
            Default True.  Draw country and coastal outlines on ``Basemap``.
        draw_latlons : Boolean
            Default True.  Draw and label lines of latitude and longitude.
        land_alpha : float
            Default 0.85.  The alpha value of the continent fill.

        """
        # Initialize
        self.mapcorners = mapcorners
        self.standard_pm = standard_pm

        self.mapcolor = mapcolor
        self.land_alpha = land_alpha
        self.coloropts = {'light': {'water': 'white',
                                    'land': '0.95'},
                          'medium': {'water': '0.625',
                                     'land': '0.775'},
                          'dark': {'water': '0.3',
                                   'land': '0.75'}}

        self.area_threshold = area_threshold
        self.resolution = {'c': '110m', 'l': '110m', 'i': '50m', 
                          'h': '10m', 'f': '10m'}.get(resolution, '110m')
        self.zborder = zborder
        self.zlandfill = zlandfill
        self.zmapbound = zmapbound
        self.zlatlon = zlatlon

        # Initialize projection
        self._set_projection(projection)

        self._set_latlonlabels(lat_labels, lon_labels)

        self.latlon_fs = latlon_fs

        self.latspacing = latlon_spacing[0]
        self.lonspacing = latlon_spacing[1]
        self.latstep = latlon_labelspacing[0]
        self.lonstep = latlon_labelspacing[1]

        self.drawstates = drawstates
        self.drawoutlines = drawoutlines
        self.draw_latlons = draw_latlons

        # Try to set label attributes
        if maplabels is not None:

            self.labels, self.labelstyle = labelfile_reader(maplabels[1])
            self.labelgroup = maplabels[0]

            # Label zorder optional, default 15 if none given.
            try:
                self.label_zorder = maplabels[2]
            except:
                self.label_zorder = 15
        else:
            self.labels = None

    def _set_latlonlabels(self, lat_labels, lon_labels):
        """
        Edit latitude, longitude labelling preferences.

        Parameters
        ----------
        lat_labels : list of strings
            The sides of the map with longitude labels.
        lon_labels : list of strings
            The sides of the map with latitude labels.

        """
        meridian_labels = [0, 0, 0, 0]
        parallel_labels = [0, 0, 0, 0]

        ind_dict = {'left': 0,
                    'right': 1,
                    'top': 2,
                    'bottom': 3}

        for la in lat_labels:
            parallel_labels[ind_dict[la]] = 1
        self.parallel_labels = parallel_labels

        for lo in lon_labels:
            meridian_labels[ind_dict[lo]] = 1
        self.meridian_labels = meridian_labels

    def _set_projection(self, projection):
        """
        Set the projection.  Defaults to 'cyl'.

        Parameters
        ----------
        projection : string
        Indicates which projection to use.  Default 'cyl'
            'cyl' : Equidistant cylindrical
            'cea' : Equal Area cylindrical
            'lcc' : Lambert Conformal Conic
            'aea' : Albers Equal Area Conic
            'ortho' : Orthographic (globe)
            'npstere' : North polar steroegraphic (conformal)
            'spstere' : South polar steroegraphic (conformal)
            'nplaea' : North polar Lambert azimuthal equal area
            'splaea' : South polar Lambert azimuthal equal area
        """
        # Map Basemap projections to Cartopy
        available_proj = {
            'cyl': ccrs.PlateCarree,
            'cea': ccrs.EqualEarth,
            'lcc': ccrs.LambertConformal,
            'aea': ccrs.AlbersEqualArea,
            'ortho': ccrs.Orthographic,
            'npstere': ccrs.NorthPolarStereo,
            'spstere': ccrs.SouthPolarStereo,
            'nplaea': lambda: ccrs.LambertAzimuthalEqualArea(central_longitude=self.standard_pm[0], central_latitude=90),
            'splaea': lambda: ccrs.LambertAzimuthalEqualArea(central_longitude=self.standard_pm[0], central_latitude=-90)
        }

        if projection in available_proj:
            if projection in ['lcc', 'aea']:
                self.proj = available_proj[projection](
                    central_longitude=self.standard_pm[0],
                    central_latitude=self.standard_pm[1],
                    standard_parallels=(self.standard_pm[2], self.standard_pm[3])
                )
            elif projection in ['npstere', 'spstere']:
                self.proj = available_proj[projection](
                    central_longitude=self.standard_pm[0]
                )
            elif projection == 'ortho':
                self.proj = available_proj[projection](
                    central_longitude=self.standard_pm[0],
                    central_latitude=self.standard_pm[1]
                )
            elif projection in ['nplaea', 'splaea']:
                self.proj = available_proj[projection]()
            else:
                self.proj = available_proj[projection]()
        else:
            self.proj = ccrs.PlateCarree()
            print('Projection not recognized, defaulting to PlateCarree.')

    def make_basemap(self, ax=None, figsize=(10, 10)):
        """Create a map using Cartopy."""
        if ax is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(1, 1, 1, projection=self.proj)
        else:
            # Check if ax has a projection by checking its type
            if not hasattr(ax, 'projection') or not isinstance(ax.projection, ccrs.Projection):
                # If ax is provided but doesn't have a projection, create new axis with projection
                fig = ax.figure
                position = ax.get_position()
                fig.delaxes(ax)
                ax = fig.add_axes(position, projection=self.proj)
        
        # Set map extent
        if hasattr(self, 'mapcorners'):
            ax.set_extent([self.mapcorners[0], self.mapcorners[2],
                          self.mapcorners[1], self.mapcorners[3]], 
                         crs=ccrs.PlateCarree())

        # Add features
        if self.mapcolor is not None:
            colors = self.coloropts[self.mapcolor]
            ax.set_facecolor(colors['water'])
            ax.add_feature(cfeature.LAND, facecolor=colors['land'],
                          alpha=self.land_alpha, zorder=self.zlandfill)
            
        if self.drawoutlines:
            ax.add_feature(cfeature.COASTLINE, zorder=self.zborder)
            ax.add_feature(cfeature.BORDERS, zorder=self.zborder)
            
        if self.drawstates:
            ax.add_feature(cfeature.STATES, zorder=self.zborder)

        # Add gridlines
        if self.draw_latlons:
            gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                            linewidth=1, color='gray', alpha=0.5,
                            linestyle='--', zorder=self.zlatlon)
            
            gl.xlocator = plt.MultipleLocator(self.lonstep)
            gl.ylocator = plt.MultipleLocator(self.latstep)
            gl.xlabel_style = {'size': self.latlon_fs}
            gl.ylabel_style = {'size': self.latlon_fs}
            
            # Set label visibility based on user preferences
            gl.top_labels = bool(self.meridian_labels[2])
            gl.bottom_labels = bool(self.meridian_labels[3])
            gl.left_labels = bool(self.parallel_labels[0])
            gl.right_labels = bool(self.parallel_labels[1])

        return ax




