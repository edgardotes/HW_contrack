#plotting functions
def plot_three_maps(bias, da_4years, da_era5_4years,
                    abs_levels=np.arange(-40, 40, 4),
                    abs_cmap='RdBu_r',
                    bias_cmap='RdBu_r',
                    titles=('Bias', '4-year mean (ICON)', '13-year mean (ERA5)')):
    """
    Plot three global maps on Robinson projection:
      1) bias (diverging, centered at 0)
      2) da_4years (absolute field, shared scale with #3)
      3) da_era5_4years (absolute field, shared scale with #2)

    Parameters
    ----------
    bias, da_4years, da_era5_4years : xarray.DataArray
        Lat/lon gridded fields in PlateCarree coordinates.
    abs_levels : 1D array-like
        Contourf levels for the absolute fields (2 & 3).
    abs_cmap : str
        Colormap for the absolute fields.
    bias_cmap : str
        Colormap for the bias field.
    titles : tuple of str
        Titles for (bias, da_4years, da_era5_4years).
    """
    proj_map = ccrs.Robinson()
    data_crs = ccrs.PlateCarree()

    fig, axes = plt.subplots(3, 1, figsize=(5, 6), dpi=100,
                             subplot_kw={'projection': proj_map})

    # ---- Absolute fields share the same normalization ----
    abs_norm = BoundaryNorm(abs_levels, ncolors=plt.get_cmap(abs_cmap).N, clip=True)

    # ---- Bias uses a centered diverging norm ----
    vmax = np.nanpercentile(np.abs(bias.values), 99)  # robust symmetric range
    if not np.isfinite(vmax) or vmax == 0:
        vmax = np.nanmax(np.abs(bias.values)) if np.isfinite(np.nanmax(np.abs(bias.values))) else 1.0
    bias_norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    # 3) da_era5_4years
    h2 = da_era5_4years.plot(ax=axes[0],
                             transform=data_crs,
                             cmap=abs_cmap,
                             levels=abs_levels,
                             norm=abs_norm,
                             cbar_kwargs={'shrink': 0.6, 'label': 'T (c)'}
                             )
    axes[0].set_title(titles[0])
    axes[0].coastlines()

    # 2) da_4years
    h1 = da_4years.plot(ax=axes[1],
                        transform=data_crs,
                        cmap=abs_cmap,
                        levels=abs_levels,
                        norm=abs_norm,
                        cbar_kwargs={'shrink': 0.6, 'label': 'T (C)'}
                        )
    axes[1].set_title(titles[1])
    axes[1].coastlines()
    
    # 1) Bias
    h0 = bias.plot(ax=axes[2],
                   transform=data_crs,
                   cmap=bias_cmap,
                   norm=bias_norm,
                   cbar_kwargs={'shrink': 0.6, 'label': 'Bias (C)'}
                   )
    axes[2].set_title(titles[2])
    axes[2].coastlines()

    # A little breathing room
    plt.tight_layout()
    return fig, (h0, h1, h2)