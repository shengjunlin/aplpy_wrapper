#! /usr/bin/env python
# #########################################################
# Author : Axel Donath (CfA)
# Email : https://gist.github.com/adonath
# Source: https://gist.github.com/adonath/c9a97d2f2d964ae7b9eb
# Description : Illustrate all the ds9 colormaps by
# $ python illustrate_ds9_cmap.py
# #########################################################
import numpy as np
import ds9_cmap

def grayify_colormap(cmap, mode='hsp'):
    """
    Return a grayscale version a the colormap.

    The grayscale conversion of the colormap is bases on perceived luminance of
    the colors. For the conversion either the `~skimage.color.rgb2gray` or a
    generic method called ``hsp`` [1]_ can be used. The code is loosely based
    on [2]_.


    Parameters
    ----------
    cmap : str or `~matplotlib.colors.Colormap`
        Colormap name or instance.
    mode : {'skimage, 'hsp'}
        Grayscale conversion method. Either ``skimage`` or ``hsp``.

    References
    ----------

    .. [1] Darel Rex Finley, "HSP Color Model - Alternative to HSV (HSB) and HSL"
       http://alienryderflex.com/hsp.html

    .. [2] Jake VanderPlas, "How Bad Is Your Colormap?"
       https://jakevdp.github.io/blog/2014/10/16/how-bad-is-your-colormap/
    """
    import matplotlib.pyplot as plt
    cmap = plt.cm.get_cmap(cmap)
    colors = cmap(np.arange(cmap.N))

    if mode == 'skimage':
        from skimage.color import rgb2gray
        luminance = rgb2gray(np.array([colors]))
        colors[:, :3] = luminance[0][:, np.newaxis]
    elif mode == 'hsp':
            RGB_weight = [0.299, 0.587, 0.114]
            luminance = np.sqrt(np.dot(colors[:, :3] ** 2, RGB_weight))
            colors[:, :3] = luminance[:, np.newaxis]
    else:
        raise ValueError('Not a valid grayscale conversion mode.')

    return cmap.from_list(cmap.name + "_grayscale", colors, cmap.N)


def illustrate_colormap(cmap, **kwargs):
    """
    Illustrate color distribution and perceived luminance of a colormap.

    Parameters
    ----------
    cmap : str or `~matplotlib.colors.Colormap`
        Colormap name or instance.
    kwargs : dicts
        Keyword arguments passed to `grayify_colormap`.
    """
    import matplotlib.pyplot as plt
    cmap = plt.cm.get_cmap(cmap)
    cmap_gray = grayify_colormap(cmap, **kwargs)
    figure = plt.figure(figsize=(6, 4))
    v = np.linspace(0, 1, 4 * cmap.N)

    # Show colormap
    show_cmap = figure.add_axes([0.1, 0.8, 0.8, 0.1])
    im = np.outer(np.ones(50), v)
    show_cmap.imshow(im, cmap=cmap, origin='lower')
    show_cmap.set_xticklabels([])
    show_cmap.set_yticklabels([])
    show_cmap.set_yticks([])
    show_cmap.set_title('RGB & Gray Luminance of colormap {0}'.format(cmap.name))

    # Show colormap gray
    show_cmap_gray = figure.add_axes([0.1, 0.72, 0.8, 0.09])
    show_cmap_gray.imshow(im, cmap=cmap_gray, origin='lower')
    show_cmap_gray.set_xticklabels([])
    show_cmap_gray.set_yticklabels([])
    show_cmap_gray.set_yticks([])

    # Plot RGB profiles
    plot_rgb = figure.add_axes([0.1, 0.1, 0.8, 0.6])
    plot_rgb.plot(v, [cmap(_)[0] for _ in v], color='r')
    plot_rgb.plot(v, [cmap(_)[1] for _ in v], color='g')
    plot_rgb.plot(v, [cmap(_)[2] for _ in v], color='b')
    plot_rgb.plot(v, [cmap_gray(_)[0] for _ in v], color='k', linestyle='--')
    plot_rgb.set_ylabel('Luminance')
    plot_rgb.set_ylim(-0.005, 1.005)

    plt.show()

illustrate_colormap('ds9grey')
illustrate_colormap('ds9a')
illustrate_colormap('ds9b')
illustrate_colormap('ds9bb')
illustrate_colormap('ds9he')
illustrate_colormap('ds9i8')
illustrate_colormap('ds9aips0')
illustrate_colormap('ds9heat')
illustrate_colormap('ds9cool')
illustrate_colormap('ds9rainbow')
illustrate_colormap('ds9rainbow_r')
