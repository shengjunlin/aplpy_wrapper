# aplpy_wrapper

This wrapper includes some tricks to fix the APLpy code,
which was tested with `APLpy == 1.1.1` and `pyregion == 1.2` in Python 2.
`APLpy >= 2` has changed a lot such that this wrapper may not fully work.

The `aplpy_wrapper` module provides two functions:
* `aplpy_plot()`: It creates and returns a `aplpy.FITSFigure` instance.
* `discrete_cmap()`: It generates a discrete color map from a given continuous color map.

## Usage:

To import this module, please do
```
>>> import sys
>>> sys.path.append([the path to the parent directory of aplpy_wrapper/])
>>> from aplpy_wrapper import *
```

One can create a matplotlib's `Figure` class instance, and plot data.
Then save the figure into a pdf/png file.

1. To plot a map, re-center it, and add a scalebar of 0.05pc at the distance of 140pc, please do
```
>>> from astropy.coordinates import SkyCoord
>>> from astropy import units as u
>>> target = SkyCoord('5:04:07.500', '+32:43:25.00', frame='fk5', unit=(u.hourangle, u.deg))
>>> fig_single = plt.figure(figsize=(11, 10))
>>> ap_fig = aplpy_plot(..., fig=fig_single,
                        recenter=target, h_deg=1./60, w_deg=1./60,
                        show_sbar=True, length_pc=0.05, dist_pc=140,
                        ...)
>>> ap_fig.save('map.pdf', dpi=300)  # fig_single.savefig() may work but there are some trouble on bbox.
```

Alternatively, one can also create an `astro_coor.Loc` instance for the target.
```
>>> from astro_coor import Coor
>>> target = Loc(RA='5:04:07.500', Dec='+32:43:25.00', name='L1512', dist_pc=140., frame='fk5')
...
>>> ap_fig = aplpy_plot(..., fig=fig_single,
                        recenter=target, h_deg=1./60, w_deg=1./60,
                        show_sbar=True, length_pc=0.05,  # 'dist_pc' will be read from 'target'
                        ...)
```

2. To plot an 1-by-2 grid of two maps, please do
```
>>> fig_grid = plt.figure(figsize=(20, 10))
>>> ap_sub1 = aplpy_plot(..., fig=fig_grid, subplot=(1, 2, 1), ...)
>>> ap_sub2 = aplpy_plot(..., fig=fig_grid, subplot=(1, 2, 2), ...)
>>> fig_grid.savefig('map_grid.pdf', dpi=300, bbox_inches='tight')
```

More examples using `aplpy_plot()`:

1. The `factor` parameter for plotting a map and contours:
```
>>> Jy2mJy = 1e3  # Converting pixels in Jy/beam to mJy/beam
```

Data is in Jy/beam but we want to display it in mJy/beam.
```
>>> ap_fig = aplpy_plot(data, factor=Jy2mJy, show_colorbar=True, ...)
```

Colorbar tickers are specified in mJy/beam.
```
>>> ap_fig.colorbar.set_ticks([0, 50, 100])
```

Drawing contours requires to load the data (in Jy/beam) again.
I can't hack `show_contour()` with the given `factor`.
Also, one may plot contours from a completely different data,
where it is not necessary to work with the same `factor`.
Therefore, the contour levels are specified
in the original units in the given data (Jy/beam).
```
>>> ap_fig.show_contour(data, colors='k', linewidths=1, smooth=None,
                        levels=np.array([0, 50, 100])/Jy2mJy)
```

2. Tick spacing
```
>>> ap_fig.ticks.set_xspacing()  # An angle in degree
```

3. More controls on the fonts with parameters:
`size`, `weight`, `family`, `style`, `variant`, `stretch`, and `frontproperties`,
where the default values are set by matplotlib or previously set values if
`set_font` has already been called. Global default values can be set by
editing the matplotlibrc file.
```
>>> ap_fig.tick_labels.set_font()  # Labels of the RA/Dec ticks
>>> ap_fig.color.set_font()  # Labels of the colorbar ticks
>>> ap_fig.color.set_axis_label_font()  # Labels of the colorbar
>>> ap_fig.scalebar.set_font()  # Label of the scalar bar
```

4. Load ds9 region files and adjust them:
* Load a file that has a ds9 'CIRCLE' region into a layer called `pol2`
```
>>> ap_fig.show_regions('POL2_region.reg', layer='pol2')
```
APLpy will saperate the region and the annotated text (if any)
into the `pol2` layer and the `pol2_txt` layer,
where the suffix `_txt` is hard-coded in APLpy.
We can therefore change them to different colors.
```
>>> ap_fig._layers['pol2'].artistlist[0].set_edgecolor('blue')
>>> ap_fig._layers['pol2'].artistlist[0].set_linwidth(2)
>>> ap_fig._layers['pol2_txt'].artistlist[0].set_color('k')
```

* Load a file that has multiple 'LINE' regions into a layer called `vec`
```
>>> ap_fig.show_regions('POL2_vectors.reg', layer='vec')
```
Use a FOR loop to change their colors.
```
>>> for seg in panel._layers['vec'].artistlist:
>>>     seg.set_color('red')
>>>     seg.set_linewidth('red')
```

## Instruction of `aplpy_plot()`
```
def aplpy_plot(data, fig=None, subplot=(1, 1, 1),
             hdu=0, dim=[0, 1], slices=[], north=False, factor=1,
             cmap='Blues', m=None, M=None, kwargs_colorscale=None,
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

```
