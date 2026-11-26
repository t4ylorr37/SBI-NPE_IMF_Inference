import astropy
import numpy as np
from astropy.io import fits
from pathlib import Path

hdul = astropy.io.fits.open("/Users/taylo/Downloads/manga-firefly-globalprop-v3_1_1-miles (1).fits")

hdul.info()

gal_info = hdul["GALAXY_INFO"].data
global_params = hdul["GLOBAL_PARAMETERS"].data

z = gal_info["REDSHIFT"]
mass = gal_info["PHOTOMETRIC_MASS"] 
age = global_params["LW_AGE_1Re"]
metals = global_params["LW_Z_1Re"]

good = np.isfinite(z) & np.isfinite(mass) & np.isfinite(age) & np.isfinite(metals) & (z <= 0.03)
print(good.sum())

young = good & (-0.58 <= age) & (age <= 0.25) # young populations
old = good & (0.25 < age) & (age <= 1.04) # intermediate to old populations
print("young:", young.sum())
print("old:", old.sum())

metal_poor = good & (-0.95 <= metals) & (metals <= -0.31)
metal_rich = good & (-0.31 < metals) & (metals <= 0.34)
print("metal rich:", metal_rich.sum())
print("metal poor:", metal_poor.sum())

lower_M = good & (7.85 <= mass) & (mass <= 9.75)
higher_M = good & (9.75 < mass) & (mass <= 11.75)
print("low mass", lower_M.sum())
print("high mass:", higher_M.sum())

old_capped_idx = np.flatnonzero(old)[:196]
old_capped = np.zeros_like(old, dtype=bool)
old_capped[old_capped_idx] = True

metal_rich_capped_idx = np.flatnonzero(metal_rich)[:192]
metal_rich_capped = np.zeros_like(metal_rich, dtype=bool)
metal_rich_capped[metal_rich_capped_idx] = True

higher_M_capped_idx = np.flatnonzero(higher_M)[:1513]
higher_M_capped = np.zeros_like(higher_M, dtype=bool)
higher_M_capped[higher_M_capped_idx] = True

good_age = young | old_capped
good_metal = metal_poor | metal_rich_capped
good_mass = lower_M | higher_M_capped

final_sample = good_age & good_metal & good_mass
print(final_sample.sum())