#! /usr/bin/env python
# -*- coding:utf-8 -*-
# #########################################################
# Author : Sheng-Jun Lin
# Email : shengjunlin@asiaa.sinica.edu.tw
# Description :
# This wrapper includes some tricks to fix the APLpy code,
# which was tested with APLpy == 1.1.1 and pyregion == 1.2 in Python 2.
# APLpy >= 2 has changed a lot such that this wrapper may not fully work.
#
# To use aplpy_plot(),
# >>> import sys
# >>> sys.path.append([PATH of this directory])
# >>> from aplpy_wrapper import *
#
# #########################################################
import aplpy
import pyregion
from matplotlib import rc, rcParams, style
import matplotlib as mpl
import matplotlib.pyplot as plt
# from astro_coor import Loc  ## For recenter()
from astropy.coordinates import SkyCoord  ## For recenter()
from astropy import units as u
# For ds9-type cmap
from ds9_cmap import *

rc('text', usetex=True)
rcParams.update({'mathtext.default': 'regular'})
style.use('classic')

def aplpy_plot(data, fig=None, subplot=(1, 1, 1),
             hdu=0, dim=[0, 1], slices=[], north=False,
             factor=1, offset=0, fix_dec_sparx=False,
             cmap='Blues', m=None, M=None, scale_mM=False,
             kwargs_colorscale=None,
             nan_color='white', alpha=1,
             # Colorbar
             show_colorbar=True,
             colorbar_label='',
             colorbar_label_fontsize=16,
             colorbar_ticklabel_fontsize=16,
             kwargs_colorbar=None,
             # Recenter
             recenter=None, w_deg=None, h_deg=None,
             # Axis labels and tickers
             RA_format='hh:mm:ss',
             Dec_format='dd:mm:ss',
             tick_color='black',
             tick_label_fontsize=16,
             hide_RA=False, hide_Dec=False,
             hide_axis_labels=True,
             # Scalar bar
             show_sbar=False, length_pc=None, length_au=None,
             dist_pc=None, kwargs_sbar=None,
             sbar_lw=4,
             sbar_color='white',
             sbar_label='',
             sbar_label_fontsize=18,
             sbar_label_fontweight='medium',
             # Beam size
             show_beam=False,
             beam=None, kwargs_beam=None,
             bm_ec='white', bm_fc='black', bm_lw=2,
             # Title
             title='', kwargs_title=None,
             **kwargs):

    '''
    Create and return a aplpy.FITSFigure instance.

    There are some examples:

    1. The "factor" parameter for plotting a map and contours:
    >>> Jy2mJy = 1e3  # Converting pixels in Jy/beam to mJy/beam

    Data is in Jy/beam but we want to display it in mJy/beam.
    >>> ap_fig = aplpy_plot(data, factor=Jy2mJy, show_colorbar=True, ...)

    Colorbar tickers are specified in mJy/beam.
    >>> ap_fig.colorbar.set_ticks([0, 50, 100])

    Drawing contours requires to load the data (in Jy/beam) again.
    I can't hack "show_contour()" with the given "factor".
    Also, one may plot contours from a completely different data,
    where it is not necessary to work with the same `factor`.
    Therefore, the contour levels are specified
    in the orginal units in the given data (Jy/beam).
    >>> ap_fig.show_contour(data, colors='k', linewidths=1, smooth=None,
                            levels=np.array([0, 50, 100])/Jy2mJy)

    2. Tick spacing
    >>> ap_fig.ticks.set_xspacing()  # An angle in degree

    3. More controls on the fonts with parameters:
    size, weight, family, style, variant, stretch, and frontproperties,
    where the default values are set by matplotlib or previously set values if
    set_font has already been called. Global default values can be set by
    editing the matplotlibrc file.
    >>> ap_fig.tick_labels.set_font()  # Labels of the RA/Dec ticks
    >>> ap_fig.color.set_font()  # Labels of the colorbar ticks
    >>> ap_fig.color.set_axis_label_font()  # Labels of the colorbar
    >>> ap_fig.scalebar.set_font()  # Label of the scalar bar

    4. Load ds9 region files and adjust them:
    * Load a file that has a ds9 'CIRCLE' region into a layer called `pol2`
    >>> ap_fig.show_regions('POL2_region.reg', layer='pol2')

    APLpy will saperate the region and the annotated text (if any)
    into the `pol2` layer and the `pol2_txt` layer,
    where the suffix `_txt` is hard-coded in APLpy.
    We can therefore change them to different colors.
    >>> ap_fig._layers['pol2'].artistlist[0].set_edgecolor('blue')
    >>> ap_fig._layers['pol2'].artistlist[0].set_linwidth(2)
    >>> ap_fig._layers['pol2_txt'].artistlist[0].set_color('k')

    * Load a file that has multiple 'LINE' regions into a layer called 'vec'
    >>> ap_fig.show_regions('POL2_vectors.reg', layer='vec')

    Use a for loop to change their colors.
    >>> for seg in panel._layers['vec'].artistlist:
    >>>     seg.set_color('red')
    >>>     seg.set_linewidth(1)

    * Load a file that has multiple 'TEXT' regions into a layer called 'labels'
    >>> for seg in panel._layers['labels_txt'].artistlist:
    >>>     print(seg.get_text())
    >>>     seg.set_color('blue')
    >>>     seg.set_fontsize('medium')
    >>>     seg.set_fontweight('bold')

    * Note: an unexpected behaviour is that 'POINT' regions are all stored
    in the layer with the suffix `_txt` while the main layer without the
    suffix becomes empty. Therefore, panel._layers[layer_name+'_txt'].artistlist
    is a list of <Line2D>, <Annotation>, <Line2D>, <Annotation>, ...,
    if the 'POINT' regions have texts.
    Otherwise, the list only has <Line2D> instances.
    In such cases, to change their colors, linewidthes, and colors:
    >>> ap_fig.show_regions('sources.reg', layer='sou')
    >>> for seg in panel._layers['sou_txt'].artistlist:
    >>>     seg.set_markeredgecolor('red')
    >>>     seg.set_markeredgewidth(1)
    >>>     seg.set_markersize(12)


    Parameters
    ----------

    data : see below
        The FITS file to open. The following data types can be passed:
             string
             astropy.io.fits.PrimaryHDU
             astropy.io.fits.ImageHDU
             astropy.wcs.WCS
             np.ndarray
             RGB image with AVM meta-data

    fig : ~matplotlib.figure.Figure, optional
        If specified, a subplot will be added to this existing
        matplotlib figure() instance, rather than a new figure
        being created from scratch.
        (Default value = None)

    subplot : tuple or list
        If specified, a subplot will be added at this position. If a tuple
        of three values, the tuple should contain the standard matplotlib
        subplot parameters, i.e. (ny, nx, subplot). If a list of four
        values, the list should contain [xmin, ymin, dx, dy] where xmin
        and ymin are the position of the bottom left corner of the
        subplot, and dx and dy are the width and height of the subplot
        respectively. These should all be given in units of the figure
        width and height. For example, [0.1, 0.1, 0.8, 0.8] will almost
        fill the entire figure, leaving a 10 percent margin on all sides.
        (Default value = (1, 1, 1))

    hdu : int
        By default, the image in the primary HDU is read in. If a
        different HDU is required, use this argument.
        (Default value = 0)

    dim : tuple or list, optional
        The index of the axes to use if the data has more than three
        dimensions.
        (Default value = [0, 1])

    slices : tuple or list, optional
        If a FITS file with more than two dimensions is specified,
        then these are the slices to extract. If all extra dimensions
        only have size 1, then this is not required.

    north : bool
        Whether to rotate the image so that the North Celestial
        Pole is up. Note that this option requires Montage to be
        installed.
        (Default value = False)

    factor : float
        A factor applied on the pixel values. This is a method hacking
        aplpy 1.1.1. (Default value = 1)

    offset : float
        An offset applied on the pixel values. This is a method hacking
        aplpy 1.1.1. (Default value = 1)

    cmap : str
        Colormap. (Default value = 'Blues')
        The ds9 type colormaps are allowed:
        'ds9b',
        'ds9cool',
        'ds9a',
        'ds9i8',
        'ds9aips0',
        'ds9rainbow',
        'ds9rainbow_r',
        'ds9he',
        'ds9heat'.

    m, M : float
        The display limits of the colorscale.
        (Default value = None, then they are 0.25th- to 99.75th-percentile)

    scale_mM : bool
        Replace `m` and `M` with the minimum and maximum of the map.
        (Default value = False)

    kwargs_colorscale : dict
        Additional arguments passed to FITSFigure.show_colorscale().
        (Default value = None)

    nan_color : str
        The background color.
        (Default value = 'w')

    alpha : float
        The transparency of the colorscale. (Default value = 1)

    show_colorbar : bool
        Show a colorbar. (Default value = True)

    colorbar_label : str
        The label of the colorbar. (Default value = '')

    colorbar_label_fontsize : int
        Equivalent to: ap_fig.colorbar.set_axis_label_font(size=...)
        (Default value = 16)

    colorbar_ticklabel_fontsize : int
        (Default value = 16)

    kwargs_colorbar : dict
        Additional arguments passed to FITSFigure.add_colorbar().
        (Default value = None)

    recenter : an astro_coor.Loc object,
               or an astropy.coordinates.SkyCoord object,
               or a 2-item list/tuple/1d-array of [RA_deg, Dec_deg].
        Note the numerical values of RA_deg and Dec_deg must be
        in the units of degree (Default value = None).

    w_deg : float
        Width in degree of the FOV if recenter is given.
        (Default value = None)

    h_deg : float
        Height in degree of the FOV if recenter is given.
        (Default value = None)

    RA_format : str
        (Default value = 'hh:mm:ss')

    Dec_format : str
        (Default value = 'dd:mm:ss')

    tick_color : str
        (Default value = 'black')

    tick_label_fontsize : int
        Set the font size of the tick labels.
        Equivalent to: ap_fig.tick_labels.set_font(size=...)
        (Default value = 16)

    hide_RA : bool
        Hide the RA ticklabels. (Default value = False)

    hide_Dec : bool
        Hide the Dec ticklabels. (Default value = False)

    hide_axis_labels : bool
        Hide the label "RA (J2000)" and "Dec (J2000)".
        (Default value = True)

    show_sbar : bool
        Show a scalar bar. (Default value = False)

    length_pc, length_au : float
        Either one should be given. (Default value = None)

    dist_pc : float
        Specify the distance in pc via this parameter.
        Otherwise, the distance is read from the given recenter.
        (Default value = None)

    kwargs_sbar : dict
        Additional arguments passed to FITSFigure.add_scalarbar().
        (Default value = None)

    sbar_lw : int
        Linewidth of the scalarbar. (Default value = 4)

    sbar_color : str
        Color of the scalarbar. (Default value = 'white')

    sbar_label : str
        Scalarbar label. (Default value = '')

    sbar_label_fontsize : int
        Font size of the scalarbar label.
        Equivalent to: ap_fig.scalebar.set_font(size=...)
        (Default value = 18)

    sbar_label_fontweight : str
        Font weight of the scalarbar label.
        Equivalent to: ap_fig.scalebar.set_font(weight=...)
        (Default value = 'medium')

    show_beam : bool
        Show the beam size. (Default value = False)

    beam : float/a 3-item list/tuple/1d-array
        If a beam is given, the beam info in the data header
        will be overwritten by this beam parameter.
        It can be a single number in arcsec for a circular beam,
        or a list of [bmaj_asec, bmin_asec, bpa_deg] for
        an elliptical beam.

    kwargs_beam : dict
        Additional arguments passed to FITSFigure.add_beam().
        (Default value = None)

    bm_ec : str
        Edge color of the beam (Default value = 'white').

    bm_fc : str
        Face color of the beam (Default value = 'black').

    bm_lw : int
        Linewidth of the beam edge (Default value = 2).

    title : str
        The title of the plot. (Default value = '')

    kwargs_title : dict
        Additional arguments passed to FITSFigure.set_title().
        (Default value = None, and then size=24 and va='bottom' will be applied)

    kwargs : dict
        Any additional arguments are passed on to matplotlib's Figure()
        class.
        For example, to set the figure size, use the
        figsize=(xsize, ysize) argument (where xsize and ysize are in
        inches). For more information on these additional arguments,
        see the *Optional keyword arguments* section in the documentation
        for `Figure <http://matplotlib.org/api/figure_api.html?#matplotlib.figure.Figure>`_

    Returns
    -------
    ap_fig : an aplpy.FITSFigure instance.

    '''

    # Make empty kwargs to be {}
    if not kwargs_colorscale:
        kwargs_colorscale = {}
    if not kwargs_colorbar:
        kwargs_colorbar = {}
    if not kwargs_sbar:
        kwargs_sbar = {}
    if not kwargs_title:
        kwargs_title = {'size': 24, 'va': 'bottom'}
    if not kwargs_beam:
        kwargs_beam = {}

    # Create the instance
    ap_fig = aplpy.FITSFigure(data,
            figure=fig, subplot=subplot,
            hdu=hdu, dimensions=dim, slices=slices, north=north,
            **kwargs)

    # Apply a factor to scale the pixels
    if factor != 1:
        ap_fig._data *= factor

    # Apply an offset to the pixels
    if offset != 0:
        ap_fig._data += offset

    # Colorscale and the background color
    if scale_mM:
        m = np.nanmin(ap_fig._data)
        M = np.nanmax(ap_fig._data)
    ap_fig.show_colorscale(cmap=cmap, vmin=m, vmax=M, **kwargs_colorscale)
    ap_fig.set_nan_color(nan_color)
    ap_fig.image.set_alpha(alpha)

    # Add a colorbar
    if show_colorbar:
        ap_fig.add_colorbar(**kwargs_colorbar)
        ap_fig.colorbar.set_font(size=colorbar_ticklabel_fontsize)
        ap_fig.colorbar.set_axis_label_text(colorbar_label)
        ap_fig.colorbar.set_axis_label_font(size=colorbar_label_fontsize)

    # Zoom in the data
    if recenter:
        # recenter = [RA_deg, Dec_deg]
        if isinstance(recenter, (list, tuple, np.ndarray)):
            ap_fig.recenter(float(recenter[0]),
                             float(recenter[1]),
                             width=w_deg, height=h_deg)
        else:
            try:
                # if recenter is an astropy.coordinates.SkyCoord instance
                ap_fig.recenter(recenter.ra.degree,
                                 recenter.dec.degree,
                                 width=w_deg, height=h_deg)
            except AttributeError:
                try:
                    # if recenter is an astr_coor.Loc instance
                    ap_fig.recenter(recenter.RA.all_d,
                                     recenter.Dec.all_d,
                                     width=w_deg, height=h_deg)
                except AttributeError:
                    raise RuntimeError('aplpy_plot: Please check if '
                            '"recenter" is an astro_coor.Loc object '
                            'or a astropy.coordinates.SkyCoord object. '
                            'In such cases, "from astro_coor import Loc" or '
                            '"from astropy.coordinates import SkyCoord" should '
                            'be executed!')

    # Format of the ticklabels and their color
    ap_fig.tick_labels.set_xformat(RA_format)
    ap_fig.tick_labels.set_yformat(Dec_format)
    ap_fig.tick_labels.set_font(size=tick_label_fontsize)
    ap_fig.ticks.set_color(tick_color)

    # Hide RA, Dec ticklabels
    if hide_RA:
        ap_fig.tick_labels.hide_x()
    else:
        ap_fig.tick_labels.show_x()
    if hide_Dec:
        ap_fig.tick_labels.hide_y()
    else:
        ap_fig.tick_labels.show_y()

    # Remove RA(J2000) and Dec(J2000)
    if hide_axis_labels:
        ap_fig.hide_xaxis_label()
        ap_fig.hide_yaxis_label()

    # Scalar bar
    if show_sbar:
        if not dist_pc:
            if ('Loc' in dir()) and isinstance(recenter, Loc):
                dist_pc = recenter.dist_pc
            else:
                raise RuntimeError('aplpy_plot: Please check if '
                        '"recenter" is an astro_coor.Loc object. '
                        'In such a case, "from astro_coor import Loc" '
                        'should be executed!')
        if length_pc:
            pc_cgs = 3.0857e+18
            au_cgs = 1.496e+13
            l_as = float(length_pc)*pc_cgs/au_cgs/float(dist_pc)
            l_deg = l_as/3600.
        elif length_au:
            l_as = float(length_au)/float(dist_pc)
            l_deg = l_as/3600.
        else:
            raise RuntimeError('aplpy_plot: Either length_pc '
                    'or length_au should be given!')
        ap_fig.add_scalebar(l_deg, **kwargs_sbar) #deg
        ap_fig.scalebar.set_linewidth(sbar_lw) #pt
        ap_fig.scalebar.set_color(sbar_color)
        ap_fig.scalebar.set_label(sbar_label)
        ap_fig.scalebar.set_font(size=sbar_label_fontsize)
        ap_fig.scalebar.set_font(weight=sbar_label_fontweight)

    if show_beam:
        if beam:
            if isinstance(beam, (int, float)):
                bmaj_asec = beam
                bmin_asec = beam
                bpa_deg = 0.
            else:
                bmaj_asec, bmin_asec, bpa_deg = beam
            ap_fig._header['BMAJ'] = float(bmaj_asec)/3600 # deg
            ap_fig._header['BMIN'] = float(bmin_asec)/3600 # deg
            ap_fig._header['BPA'] = float(bpa_deg) # deg
        ap_fig.add_beam(**kwargs_beam)
        ap_fig.beam.set_edgecolor(bm_ec)
        ap_fig.beam.set_facecolor(bm_fc)
        ap_fig.beam.set_linewidth(bm_lw)

    if title:
        # APLpy seems to reset the rcParams so I update it again
        rcParams.update({'mathtext.default': 'regular'})
        ap_fig.set_title(title, **kwargs_title)

    return ap_fig

def discrete_cmap(bounds, cmap_template=plt.cm.jet, first_grey=False):
    """Generate a discrete color map.

    # x and y are the data, while c and s are their colors and sizes.
    # to discrete the given camp into N bins
    bounds = np.linspace(vmin, vmax, N+1)
    cmap, norm = discrete_cmap(bounds)
    # make the scatter
    scat = ax.scatter(x, y, c, s,
                      cmap=cmap, norm=norm)
    # create a second axes for the colorbar
    ax2 = fig.add_axes([0.95, 0.1, 0.03, 0.8])
    cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm,
        spacing='proportional', ticks=bounds, boundaries=bounds, format='%1i')

    Ref: https://stackoverflow.com/questions/14777066/matplotlib-discrete-colorbar

    Parameters
    ----------
    bounds : list
        The discrete points in the colorbar.

    cmap_tamplate : plt.cm object
        The camp for discreting.

    first_grey : bool
        Force the first color entry to be grey. (Default value = False)

    Returns
    -------
    (discrete_cmap, norm)

    discrete_cmap : matplotlib.pyplot.cm object

    norm : matplotlib.colors.Normalize object

    """
    cmap = cmap_template  # define the colormap
    # extract all colors from the .jet map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    if first_grey:
        # force the first color entry to be grey
        cmaplist[0] = (.5, .5, .5, 1.0)

    # create the new map
    discrete_cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'Custom cmap', cmaplist, cmap.N)
    norm = mpl.colors.BoundaryNorm(bounds, discrete_cmap.N)
    return (discrete_cmap, norm)

