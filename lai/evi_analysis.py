#!/usr/bin/env python3
"""
EVI Comparison Analysis Module
==============================

This module calculates EVI (Enhanced Vegetation Index) for both Sentinel-2 and UAV data
and creates a comprehensive comparison plot with error margins, lab measurements, and
field measurements.

Key Features:
- Calculates EVI for both Sentinel-2 and UAV platforms
- Interpolates EVI values to field measurement dates
- Creates plots with error margins around EVI lines
- Includes lab and field measurements
- Clean plot design with no background/borders
- Custom date formatting (DD.MM.YY)

Author: Kai; created with AI support
Date: 2025
"""

import numpy as np
import pandas as pd
import rasterio
from rasterio.mask import mask
from rasterio.features import geometry_mask
import geopandas as gpd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import interp1d
from scipy.stats import pearsonr, spearmanr
from scipy.interpolate import CubicSpline
import os
import glob
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Field measurements data
FIELD_MEASUREMENTS = {
    '2025-04-08': 1.68,
    '2025-04-15': 2.98,
    '2025-04-22': 4.77,
    '2025-04-29': 3.69,
    '2025-05-06': 4.70,
    '2025-05-13': 3.48,
    '2025-05-20': 4.41,
    '2025-05-27': 5.4,
    '2025-06-03': 5.69,
    '2025-06-10': 5.81,
    '2025-06-17': 4.9,
    '2025-06-24': 3.98,
}

# LAB measurements from the field   
LAB_MEASUREMENTS = {
    '2025-04-08': 0.906904762,
    '2025-04-15': 1.474109524,
    '2025-04-22': 3.066285714,
    '2025-04-29': 2.73147619,
    '2025-05-06': 3.344857143,
    '2025-05-13': 2.663952381,
    '2025-05-20': 1.62452381,
    '2025-05-27': 1.744952381,
    '2025-06-03': 1.806380952,
    '2025-06-10': 1.0979,
    '2025-06-24': 0,
}

class EVIAnalyzer:
    """
    EVI analysis for both Sentinel-2 and UAV platforms
    """
    
    def __init__(self, sentinel_input_dir, uav_input_dir, shapes_file, output_dir):
        """
        Initialize the EVI analyzer
        
        Parameters:
        -----------
        sentinel_input_dir : str
            Directory containing Sentinel-2 red-edge images
        uav_input_dir : str
            Directory containing UAV multispectral images
        shapes_file : str
            Path to GeoJSON file with field boundaries
        output_dir : str
            Output directory for results
        """
        self.sentinel_input_dir = sentinel_input_dir
        self.uav_input_dir = uav_input_dir
        self.shapes_file = shapes_file
        self.output_dir = output_dir
        
        # Create output directories
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(os.path.join(output_dir, 'evi_data'), exist_ok=True)
        os.makedirs(os.path.join(output_dir, 'lai_data'), exist_ok=True)
        os.makedirs(os.path.join(output_dir, 'plots'), exist_ok=True)
        os.makedirs(os.path.join(output_dir, 'cache'), exist_ok=True)
        
        # Load field boundaries
        self.shapes = self.load_shapes()
        
        print(f"Initialized EVI Analyzer")
        print(f"Sentinel-2 input: {sentinel_input_dir}")
        print(f"UAV input: {uav_input_dir}")
        print(f"Output directory: {output_dir}")
        print(f"Field boundaries: {len(self.shapes)} polygons loaded")
    
    def load_shapes(self):
        """Load field boundary polygons"""
        try:
            shapes = gpd.read_file(self.shapes_file)
            print(f"Loaded {len(shapes)} field boundary polygons")
            return shapes
        except Exception as e:
            print(f"Error loading shapes: {e}")
            return None
    
    def calculate_sentinel_evi(self, image_path):
        """
        Calculate EVI for Sentinel-2 image
        
        Parameters:
        -----------
        image_path : str
            Path to Sentinel-2 red-edge image
            
        Returns:
        --------
        dict
            Dictionary containing EVI values and statistics
        """
        print(f"  Calculating EVI for Sentinel-2: {os.path.basename(image_path)}...")
        
        with rasterio.open(image_path) as src:
            # Read bands (B4, B5, B6, B7, B8A, SCL)
            bands = src.read()
            
            # Handle NoData values for Sentinel-2 (common values: 0, -9999, -32768)
            bands = bands.astype(np.float32)
            nodata_values = [src.nodata, 0, -9999, -32768] if src.nodata is not None else [0, -9999, -32768]
            for nodata_val in nodata_values:
                if nodata_val is not None:
                    bands = np.where(bands == nodata_val, np.nan, bands)
            
            # Extract bands for EVI calculation
            red = bands[0]   # B4 - Red (665 nm)
            b8 = bands[4]    # B8A - NIR (865nm)
            
            # For Sentinel-2, handle blue band properly
            # Check if we have enough bands for blue (B2 would be band 6 if it exists)
            if bands.shape[0] >= 7:
                blue = bands[6]  # B2 - Blue (490 nm)
            else:
                # If no blue band, use a small constant value to avoid division issues
                blue = np.full_like(red, 0.1)  # Small constant value
            
            # Calculate EVI
            evi_numerator = (b8 - red)
            evi_denominator = (b8 + 6 * red - 7.5 * blue) + 1
            evi = np.where(evi_denominator != 0, 2.5 * evi_numerator / evi_denominator, np.nan)
            
            # Handle invalid values
            evi = np.where(np.isfinite(evi), evi, np.nan)
            
            return {
                'EVI': evi,
                'profile': src.profile.copy()
            }
    
    def calculate_uav_evi(self, image_path):
        """
        Calculate EVI for UAV image
        
        Parameters:
        -----------
        image_path : str
            Path to UAV multispectral image
            
        Returns:
        --------
        dict
            Dictionary containing EVI values and statistics
        """
        print(f"  Calculating EVI for UAV: {os.path.basename(image_path)}...")
        
        # Clip to shapes first
        clipped_bands, clipped_meta = self._clip_to_shapes(image_path)
        
        # After clipping: order is [Blue, Red, RedEdge, NIR]
        blue = clipped_bands[0].astype(np.float32)
        red = clipped_bands[1].astype(np.float32)
        nir = clipped_bands[3].astype(np.float32)  # NIR band
        
        print(f"    UAV band shapes: Blue={blue.shape}, Red={red.shape}, NIR={nir.shape}")
        print(f"    UAV band ranges: Blue=[{np.nanmin(blue):.6f}, {np.nanmax(blue):.6f}], Red=[{np.nanmin(red):.6f}, {np.nanmax(red):.6f}], NIR=[{np.nanmin(nir):.6f}, {np.nanmax(nir):.6f}]")
        print(f"    UAV valid pixels per band: Blue={np.sum(np.isfinite(blue))}, Red={np.sum(np.isfinite(red))}, NIR={np.sum(np.isfinite(nir))}")
        
        # Calculate EVI - only where all bands have valid values
        valid_mask = np.isfinite(blue) & np.isfinite(red) & np.isfinite(nir)
        
        # Initialize EVI array with NaN
        evi = np.full_like(blue, np.nan)
        
        # Calculate EVI only for valid pixels
        if np.sum(valid_mask) > 0:
            blue_valid = blue[valid_mask]
            red_valid = red[valid_mask]
            nir_valid = nir[valid_mask]
            
            evi_denominator = (nir_valid + 6 * red_valid - 7.5 * blue_valid) + 1
            evi_valid = np.where(evi_denominator != 0, 2.5 * ((nir_valid - red_valid) / evi_denominator), np.nan)
            
            # Put valid EVI values back into the full array
            evi[valid_mask] = evi_valid
        
        print(f"    UAV EVI range: [{np.nanmin(evi):.6f}, {np.nanmax(evi):.6f}], valid pixels: {np.sum(np.isfinite(evi))}")
        
        return {
            'EVI': evi,
            'profile': clipped_meta
        }
    
    def _clip_to_shapes(self, image_path):
        """
        Clip image to field boundaries
        
        Parameters:
        -----------
        image_path : str
            Path to image file
            
        Returns:
        --------
        tuple
            (clipped_bands, metadata)
        """
        with rasterio.open(image_path) as src:
            # Read all bands
            bands = src.read()
            print(f"    Original UAV image shape: {bands.shape}")
            print(f"    Original data type: {bands.dtype}")
            print(f"    NoData value: {src.nodata}")
            
            # Handle NoData values - common UAV NoData values are -32767, -9999, 0
            nodata_values = [src.nodata, -32767, -9999] if src.nodata is not None else [-32767, -9999]
            
            # Convert to float and handle NoData
            bands = bands.astype(np.float32)
            for nodata_val in nodata_values:
                if nodata_val is not None:
                    bands = np.where(bands == nodata_val, np.nan, bands)
            
            print(f"    After NoData handling - value ranges: min={np.nanmin(bands):.6f}, max={np.nanmax(bands):.6f}")
            
            # Create mask for all polygons
            mask_geometry = [polygon.geometry for polygon in self.shapes.itertuples()]
            
            if mask_geometry:
                try:
                    # Clip the image
                    clipped_bands, clipped_transform = mask(src, mask_geometry, crop=True, nodata=np.nan)
                    print(f"    Clipped UAV image shape: {clipped_bands.shape}")
                    
                    # Additional NoData handling after clipping
                    clipped_bands = clipped_bands.astype(np.float32)
                    for nodata_val in nodata_values:
                        if nodata_val is not None:
                            clipped_bands = np.where(clipped_bands == nodata_val, np.nan, clipped_bands)
                    
                    print(f"    After clipping - value ranges: min={np.nanmin(clipped_bands):.6f}, max={np.nanmax(clipped_bands):.6f}")
                    
                    clipped_meta = src.meta.copy()
                    clipped_meta.update({
                        'height': clipped_bands.shape[1],
                        'width': clipped_bands.shape[2],
                        'transform': clipped_transform,
                        'dtype': 'float32',
                        'nodata': np.nan
                    })
                    return clipped_bands, clipped_meta
                except Exception as e:
                    print(f"    Warning: Could not clip image: {e}")
                    return bands, src.meta
            else:
                print(f"    Warning: No mask geometry found, using original image")
                return bands, src.meta
    
    def extract_zonal_statistics(self, image_path, evi_data, platform):
        """
        Extract zonal statistics for EVI
        
        Parameters:
        -----------
        image_path : str
            Path to the image file
        evi_data : dict
            Dictionary containing EVI values
        platform : str
            Platform name ('Sentinel-2' or 'UAV')
            
        Returns:
        --------
        list
            List of dictionaries with zonal statistics
        """
        print("  Extracting zonal statistics...")
        
        # Extract date from filename
        filename = os.path.basename(image_path)
        if platform == 'Sentinel-2':
            date_str = filename.split('_')[2]  # Extract date from filename
        else:  # UAV
            date_str = filename[:8]  # Extract YYYYMMDD
            
        zonal_results = []
        
        # Get EVI values and profile
        evi_values = evi_data['EVI']
        profile = evi_data['profile']
        
        # For UAV images, process each polygon individually like Sentinel-2
        if platform == 'UAV':
            # Use the same polygon-by-polygon approach as Sentinel-2
            # We need to work with the original unclipped image to get proper polygon masks
            original_image_path = image_path.replace('_clipped', '')  # In case it was clipped
            
            # Read the original UAV image to get proper coordinate system
            with rasterio.open(image_path) as src:
                # Get the profile from the clipped data
                clipped_profile = evi_data['profile']
                
                # Process each polygon individually
                for idx, polygon in self.shapes.iterrows():
                    try:
                        # Build a boolean mask for the polygon footprint using clipped image coordinates
                        poly_mask = geometry_mask([polygon.geometry], 
                                                  transform=clipped_profile['transform'],
                                                  invert=True, 
                                                  out_shape=(clipped_profile['height'], clipped_profile['width']))
                        
                        # Extract EVI statistics within the polygon
                        valid_values = evi_values[poly_mask]
                        valid_values = valid_values[np.isfinite(valid_values)]
                        
                        if len(valid_values) > 0:
                            polygon_stats = {
                                'date': date_str,
                                'platform': platform,
                                'polygon_id': idx,
                                'evi_mean': np.mean(valid_values),
                                'evi_std': np.std(valid_values),
                                'evi_min': np.min(valid_values),
                                'evi_max': np.max(valid_values),
                                'evi_count': len(valid_values)
                            }
                            print(f"    UAV polygon {idx}: mean={polygon_stats['evi_mean']:.3f}, std={polygon_stats['evi_std']:.3f}, count={polygon_stats['evi_count']}")
                        else:
                            polygon_stats = {
                                'date': date_str,
                                'platform': platform,
                                'polygon_id': idx,
                                'evi_mean': np.nan,
                                'evi_std': np.nan,
                                'evi_min': np.nan,
                                'evi_max': np.nan,
                                'evi_count': 0
                            }
                            print(f"    UAV polygon {idx}: No valid values")
                        
                        zonal_results.append(polygon_stats)
                        
                    except Exception as e:
                        print(f"    Error processing UAV polygon {idx}: {e}")
                        # Add a NaN record to maintain consistency
                        polygon_stats = {
                            'date': date_str,
                            'platform': platform,
                            'polygon_id': idx,
                            'evi_mean': np.nan,
                            'evi_std': np.nan,
                            'evi_min': np.nan,
                            'evi_max': np.nan,
                            'evi_count': 0
                        }
                        zonal_results.append(polygon_stats)
                        continue
        
        else:  # Sentinel-2 - use original approach with polygon masks
            with rasterio.open(image_path) as src:
                for idx, polygon in self.shapes.iterrows():
                    try:
                        # Build a boolean mask for the polygon footprint
                        poly_mask = geometry_mask([polygon.geometry], transform=src.transform,
                                                  invert=True, out_shape=(src.height, src.width))
                        
                        # Extract EVI statistics within the polygon
                        valid_values = evi_values[poly_mask]
                        valid_values = valid_values[np.isfinite(valid_values)]
                        
                        if len(valid_values) > 0:
                            polygon_stats = {
                                'date': date_str,
                                'platform': platform,
                                'polygon_id': idx,
                                'evi_mean': np.mean(valid_values),
                                'evi_std': np.std(valid_values),
                                'evi_min': np.min(valid_values),
                                'evi_max': np.max(valid_values),
                                'evi_count': len(valid_values)
                            }
                        else:
                            polygon_stats = {
                                'date': date_str,
                                'platform': platform,
                                'polygon_id': idx,
                                'evi_mean': np.nan,
                                'evi_std': np.nan,
                                'evi_min': np.nan,
                                'evi_max': np.nan,
                                'evi_count': 0
                            }
                        
                        zonal_results.append(polygon_stats)
                        
                    except Exception as e:
                        print(f"    Error processing polygon {idx}: {e}")
                        continue
        
        return zonal_results
    
    def save_platform_data(self, platform, zonal_results):
        """
        Save platform-specific EVI and LAI data to cache
        
        Parameters:
        -----------
        platform : str
            Platform name ('sentinel' or 'uav')
        zonal_results : list
            List of zonal analysis results
        """
        if not zonal_results:
            return
            
        # Convert to DataFrame
        df = pd.DataFrame(zonal_results)
        
        # Save to cache
        cache_path = os.path.join(self.output_dir, 'cache', f'{platform}_evi_lai_data.csv')
        df.to_csv(cache_path, index=False)
        print(f"  ✓ Cached {platform} data: {cache_path}")
        
        # Also save to individual folders for user access
        evi_path = os.path.join(self.output_dir, 'evi_data', f'{platform}_evi_values.csv')
        lai_path = os.path.join(self.output_dir, 'lai_data', f'{platform}_lai_values.csv')
        
        # Extract EVI and LAI specific columns
        evi_cols = ['date', 'platform', 'polygon_id'] + [col for col in df.columns if 'evi' in col.lower()]
        lai_cols = ['date', 'platform', 'polygon_id'] + [col for col in df.columns if 'lai' in col.lower()]
        
        if len(evi_cols) > 3:  # More than just date, platform, polygon_id
            df[evi_cols].to_csv(evi_path, index=False)
            print(f"  ✓ Saved {platform} EVI values: {evi_path}")
            
        if len(lai_cols) > 3:  # More than just date, platform, polygon_id  
            df[lai_cols].to_csv(lai_path, index=False)
            print(f"  ✓ Saved {platform} LAI values: {lai_path}")
    
    def load_platform_data(self, platform):
        """
        Load platform-specific EVI and LAI data from cache
        
        Parameters:
        -----------
        platform : str
            Platform name ('sentinel' or 'uav')
            
        Returns:
        --------
        list or None
            List of zonal analysis results or None if not cached
        """
        cache_path = os.path.join(self.output_dir, 'cache', f'{platform}_evi_lai_data.csv')
        
        if os.path.exists(cache_path):
            try:
                df = pd.read_csv(cache_path)
                print(f"  ✓ Loaded cached {platform} data: {len(df)} records")
                return df.to_dict('records')
            except Exception as e:
                print(f"  ⚠️ Error loading cached {platform} data: {e}")
                return None
        else:
            print(f"  ℹ️ No cached {platform} data found")
            return None
    
    def clear_cache(self, platform=None):
        """
        Clear cached data for specified platform or all platforms
        
        Parameters:
        -----------
        platform : str or None
            Platform name ('sentinel', 'uav') or None for all platforms
        """
        cache_dir = os.path.join(self.output_dir, 'cache')
        
        if platform:
            cache_file = os.path.join(cache_dir, f'{platform}_evi_lai_data.csv')
            if os.path.exists(cache_file):
                os.remove(cache_file)
                print(f"✓ Cleared {platform} cache")
        else:
            # Clear all cache files
            for file in glob.glob(os.path.join(cache_dir, '*_evi_lai_data.csv')):
                os.remove(file)
                print(f"✓ Cleared cache: {os.path.basename(file)}")
    
    def process_sentinel_images(self):
        """
        Process all Sentinel-2 images and calculate EVI (with caching)
        """
        print("Processing Sentinel-2 images for EVI calculation...")
        
        # Try to load cached data first
        cached_results = self.load_platform_data('sentinel')
        if cached_results:
            print("Using cached Sentinel-2 data, skipping processing.")
            return cached_results
        
        # Find all image files
        image_pattern = os.path.join(self.sentinel_input_dir, 'S2_RedEdge_*.tif')
        image_files = glob.glob(image_pattern)
        image_files.sort()
        
        if not image_files:
            print(f"No Sentinel-2 red-edge images found in {self.sentinel_input_dir}")
            return None
        
        print(f"Found {len(image_files)} Sentinel-2 images to process")
        
        all_zonal_results = []
        
        for i, image_path in enumerate(image_files, 1):
            print(f"\n[{i}/{len(image_files)}] Processing: {os.path.basename(image_path)}")
            
            try:
                # Calculate EVI
                evi_data = self.calculate_sentinel_evi(image_path)
                
                # Extract zonal statistics
                zonal_results = self.extract_zonal_statistics(image_path, evi_data, 'Sentinel-2')
                all_zonal_results.extend(zonal_results)
                
                print(f"  ✓ Successfully processed: {os.path.basename(image_path)}")
                
            except Exception as e:
                print(f"  ✗ Error processing {image_path}: {e}")
                continue
        
        # Save results to cache
        self.save_platform_data('sentinel', all_zonal_results)
        
        return all_zonal_results
    
    def process_uav_images(self):
        """
        Process all UAV images and calculate EVI (with caching)
        """
        print("Processing UAV images for EVI calculation...")
        
        # Try to load cached data first
        cached_results = self.load_platform_data('uav')
        if cached_results:
            print("Using cached UAV data, skipping processing.")
            return cached_results
        
        # Find all UAV image files
        image_pattern = os.path.join(self.uav_input_dir, '*_altum_ortho_norm_5cm.tif')
        image_files = glob.glob(image_pattern)
        image_files.sort()
        
        if not image_files:
            print(f"No UAV multispectral images found in {self.uav_input_dir}")
            return None
        
        print(f"Found {len(image_files)} UAV images to process")
        
        all_zonal_results = []
        
        for i, image_path in enumerate(image_files, 1):
            print(f"\n[{i}/{len(image_files)}] Processing: {os.path.basename(image_path)}")
            
            try:
                # Calculate EVI
                evi_data = self.calculate_uav_evi(image_path)
                
                # Extract zonal statistics
                zonal_results = self.extract_zonal_statistics(image_path, evi_data, 'UAV')
                all_zonal_results.extend(zonal_results)
                
                print(f"  ✓ Successfully processed: {os.path.basename(image_path)}")
                
            except Exception as e:
                print(f"  ✗ Error processing {image_path}: {e}")
                continue
        
        # Save results to cache
        self.save_platform_data('uav', all_zonal_results)
        
        return all_zonal_results
    
    def interpolate_to_field_dates(self, sentinel_results, uav_results, interpolation_method='spline'):
        """
        Interpolate EVI values to daily intervals, convert to LAI, and merge with field measurements
        
        Parameters:
        -----------
        sentinel_results : list
            List of Sentinel-2 zonal results
        uav_results : list
            List of UAV zonal results
        interpolation_method : str
            'linear' or 'spline' interpolation method
            
        Returns:
        --------
        pandas.DataFrame
            DataFrame with interpolated LAI values converted from EVI
        """
        print(f"\nInterpolating EVI data to daily intervals using {interpolation_method} method and converting to LAI...")
        
        # EVI to LAI conversion coefficients (from LAI_COEFFICIENTS)
        EVI_TO_LAI_A = 3.618
        EVI_TO_LAI_B = -0.118
        
        # Convert results to DataFrames
        sentinel_df = pd.DataFrame(sentinel_results) if sentinel_results else pd.DataFrame()
        uav_df = pd.DataFrame(uav_results) if uav_results else pd.DataFrame()
        
        # Create 4-day interval dates from the earliest to latest available dates
        all_dates = []
        
        # Add Sentinel-2 dates
        if not sentinel_df.empty:
            sentinel_dates = pd.to_datetime(sentinel_df['date'], format='%Y-%m-%d', errors='coerce')
            all_dates.extend(sentinel_dates.dropna())
        
        # Add UAV dates (format: YYYYMMDD)
        if not uav_df.empty:
            uav_dates = pd.to_datetime(uav_df['date'], format='%Y%m%d', errors='coerce')
            all_dates.extend(uav_dates.dropna())
        
        # Add field measurement dates
        field_dates = [pd.to_datetime(date_str) for date_str in FIELD_MEASUREMENTS.keys()]
        all_dates.extend(field_dates)
        
        if not all_dates:
            print("No dates available for interpolation")
            return pd.DataFrame()
        
        # Create 4-day interval dates
        min_date = min(all_dates)
        max_date = max(all_dates)
        
        # Generate daily intervals using spline interpolation
        interpolation_dates = pd.date_range(start=min_date, end=max_date, freq='D')
        
        print(f"  Interpolation range: {min_date} to {max_date}")
        print(f"  Generated {len(interpolation_dates)} daily interpolation dates")
        
        # Debug: Check if 2025-06-17 is in interpolation dates
        target_date = pd.to_datetime('2025-06-17')
        if target_date in interpolation_dates:
            print(f"  ✅ 2025-06-17 is in interpolation grid")
        else:
            print(f"  ❌ 2025-06-17 is NOT in interpolation grid")
        
        # Process Sentinel-2 data
        sentinel_interpolated = None
        if not sentinel_df.empty:
            print(f"Processing Sentinel-2 data: {len(sentinel_df)} records")
            # Calculate mean across polygons for each date
            sentinel_data = sentinel_df.groupby('date')['evi_mean'].agg(['mean', 'std']).reset_index()
            sentinel_data.columns = ['date', 'sentinel_evi_mean', 'sentinel_evi_std']
            
            # Convert EVI to LAI using formula: LAI = 3.618 × EVI - 0.118
            sentinel_data['sentinel_lai_mean'] = EVI_TO_LAI_A * sentinel_data['sentinel_evi_mean'] + EVI_TO_LAI_B
            sentinel_data['sentinel_lai_std'] = EVI_TO_LAI_A * sentinel_data['sentinel_evi_std']  # Scale std by coefficient
            
            # Convert dates
            sentinel_data['date'] = pd.to_datetime(sentinel_data['date'], format='%Y-%m-%d', errors='coerce')
            sentinel_dates_numeric = sentinel_data['date'].astype(np.int64) // 10**9
            
            # Interpolate to 4-day intervals across the full range
            if len(sentinel_dates_numeric) > 0:
                interp_dates_numeric = interpolation_dates.astype(np.int64) // 10**9
                
                # Interpolate across the full interpolation range (not just Sentinel range)
                # This allows Sentinel data to be extrapolated to UAV dates and vice versa
                if len(interp_dates_numeric) > 0:
                    # Interpolate LAI mean and std using specified method
                    if interpolation_method == 'spline' and len(sentinel_dates_numeric) >= 4:
                        interp_mean = interp1d(sentinel_dates_numeric, sentinel_data['sentinel_lai_mean'], 
                                              kind='cubic', fill_value='extrapolate', bounds_error=False)
                        interp_std = interp1d(sentinel_dates_numeric, sentinel_data['sentinel_lai_std'], 
                                             kind='cubic', fill_value='extrapolate', bounds_error=False)
                        method_used = 'cubic spline'
                    else:  # Use linear interpolation
                        interp_mean = interp1d(sentinel_dates_numeric, sentinel_data['sentinel_lai_mean'], 
                                              kind='linear', fill_value='extrapolate', bounds_error=False)
                        interp_std = interp1d(sentinel_dates_numeric, sentinel_data['sentinel_lai_std'], 
                                             kind='linear', fill_value='extrapolate', bounds_error=False)
                        method_used = 'linear'
                    
                    sentinel_interpolated = pd.DataFrame({
                        'date': interpolation_dates,
                        'sentinel_lai_mean': interp_mean(interp_dates_numeric),
                        'sentinel_lai_std': interp_std(interp_dates_numeric)
                    })
                    print(f"  Interpolated Sentinel-2 to {len(sentinel_interpolated)} daily points using {method_used}")
        
        # Process UAV data
        uav_interpolated = None
        if not uav_df.empty:
            print(f"Processing UAV data: {len(uav_df)} records")
            
            # Debug UAV EVI data
            print(f"  UAV DataFrame columns: {uav_df.columns.tolist()}")
            print(f"  UAV polygon distribution per date:")
            polygon_counts = uav_df.groupby('date')['polygon_id'].count()
            print(f"    Polygons per date: {polygon_counts.unique()}")
            print(f"    Expected: 4 polygons per date")
            print(f"  UAV EVI data preview:")
            print(uav_df[['date', 'polygon_id', 'evi_mean', 'evi_std']].head(8))
            print(f"  UAV EVI stats: min={uav_df['evi_mean'].min():.3f}, max={uav_df['evi_mean'].max():.3f}, mean={uav_df['evi_mean'].mean():.3f}")
            
            # Calculate mean across polygons for each date (UAV now has 4 polygons per date like Sentinel)
            uav_data = uav_df.groupby('date')['evi_mean'].agg(['mean', 'std']).reset_index()
            uav_data.columns = ['date', 'uav_evi_mean', 'uav_evi_std']
            
            print(f"  UAV aggregated EVI data:")
            print(uav_data.head())
            
            # Convert EVI to LAI using formula: LAI = 3.618 × EVI - 0.118
            uav_data['uav_lai_mean'] = EVI_TO_LAI_A * uav_data['uav_evi_mean'] + EVI_TO_LAI_B
            uav_data['uav_lai_std'] = EVI_TO_LAI_A * uav_data['uav_evi_std']  # Scale std by coefficient
            
            print(f"  UAV LAI conversion results:")
            print(uav_data[['date', 'uav_evi_mean', 'uav_lai_mean']].head())
            
            # Convert dates (UAV format: YYYYMMDD)
            uav_data['date'] = pd.to_datetime(uav_data['date'], format='%Y%m%d', errors='coerce')
            print(f"  UAV date range: {uav_data['date'].min()} to {uav_data['date'].max()}")
            
            # Debug: Check for specific date 2025-06-17
            target_date = pd.to_datetime('2025-06-17')
            if target_date in uav_data['date'].values:
                target_row = uav_data[uav_data['date'] == target_date]
                print(f"  ✅ Found UAV data for 2025-06-17: LAI={target_row['uav_lai_mean'].iloc[0]:.3f}")
            else:
                print(f"  ❌ UAV data for 2025-06-17 not found in processed data")
            
            uav_dates_numeric = uav_data['date'].astype(np.int64) // 10**9
            
            # Interpolate to 4-day intervals across the full range
            if len(uav_dates_numeric) > 0:
                interp_dates_numeric = interpolation_dates.astype(np.int64) // 10**9
                
                # Interpolate across the full interpolation range (not just UAV range)
                # This allows UAV data to be extrapolated to Sentinel dates and vice versa
                if len(interp_dates_numeric) > 0:
                    # Interpolate LAI mean and std using specified method
                    if interpolation_method == 'spline' and len(uav_dates_numeric) >= 4:
                        interp_mean = interp1d(uav_dates_numeric, uav_data['uav_lai_mean'], 
                                              kind='cubic', fill_value='extrapolate', bounds_error=False)
                        interp_std = interp1d(uav_dates_numeric, uav_data['uav_lai_std'], 
                                             kind='cubic', fill_value='extrapolate', bounds_error=False)
                        method_used = 'cubic spline'
                    else:  # Use linear interpolation
                        interp_mean = interp1d(uav_dates_numeric, uav_data['uav_lai_mean'], 
                                              kind='linear', fill_value='extrapolate', bounds_error=False)
                        interp_std = interp1d(uav_dates_numeric, uav_data['uav_lai_std'], 
                                             kind='linear', fill_value='extrapolate', bounds_error=False)
                        method_used = 'linear'
                    
                    uav_interpolated = pd.DataFrame({
                        'date': interpolation_dates,
                        'uav_lai_mean': interp_mean(interp_dates_numeric),
                        'uav_lai_std': interp_std(interp_dates_numeric)
                    })
                    print(f"  Interpolated UAV to {len(uav_interpolated)} daily points using {method_used}")
                    
                    # Debug: Check if 2025-06-17 is in interpolated data
                    target_date = pd.to_datetime('2025-06-17')
                    if target_date in uav_interpolated['date'].values:
                        target_row = uav_interpolated[uav_interpolated['date'] == target_date]
                        print(f"  ✅ 2025-06-17 found in interpolated UAV data: LAI={target_row['uav_lai_mean'].iloc[0]:.3f}")
                    else:
                        print(f"  ❌ 2025-06-17 NOT found in interpolated UAV data")
        
        # Create base results with field measurements
        field_dates = [pd.to_datetime(date_str) for date_str in FIELD_MEASUREMENTS.keys()]
        results = pd.DataFrame({
            'date': field_dates,
            'field_lai': list(FIELD_MEASUREMENTS.values()),
            'lab_lai': [LAB_MEASUREMENTS.get(date_str, np.nan) for date_str in FIELD_MEASUREMENTS.keys()]
        })
        
        # Normalize all dates to midnight before merging to avoid timestamp precision issues
        results['date'] = pd.to_datetime(results['date']).dt.floor('D')
        
        # Combine all unique dates from field measurements and interpolated data
        all_valid_dates = set(results['date'])
        
        if sentinel_interpolated is not None:
            sentinel_interpolated['date'] = pd.to_datetime(sentinel_interpolated['date']).dt.floor('D')
            all_valid_dates.update(sentinel_interpolated['date'])
            results = pd.merge(results, sentinel_interpolated, on='date', how='outer')
        
        if uav_interpolated is not None:
            uav_interpolated['date'] = pd.to_datetime(uav_interpolated['date']).dt.floor('D')
            all_valid_dates.update(uav_interpolated['date'])
            results = pd.merge(results, uav_interpolated, on='date', how='outer')
        
        # Final date normalization after all merges
        results['date'] = pd.to_datetime(results['date']).dt.floor('D')
        
        # Remove any rows with invalid dates (NaT values that become 1970-01-01)
        print(f"  Before filtering: {len(results)} rows")
        print(f"  Date range before filtering: {results['date'].min()} to {results['date'].max()}")
        
        # Check for problematic dates
        problematic_dates = results[results['date'].isna() | (results['date'] < pd.Timestamp('2000-01-01'))]
        if len(problematic_dates) > 0:
            print(f"  Found {len(problematic_dates)} problematic dates:")
            print(problematic_dates['date'].head())
        
        results = results[results['date'].notna()]
        results = results[results['date'] >= pd.Timestamp('2000-01-01')]  # Remove any dates before 2000
        
        print(f"  After filtering: {len(results)} rows")
        print(f"  Date range after filtering: {results['date'].min()} to {results['date'].max()}")
        
        # Sort by date
        results = results.sort_values('date').reset_index(drop=True)
        
        print(f"Successfully interpolated to {len(results)} daily dates using {interpolation_method} method")
        
        return results
    
    def create_evi_comparison_plot(self, evi_data, interpolation_method='spline'):
        """
        Create LAI comparison plot with 95% prediction intervals, lab and field measurements
        
        Parameters:
        -----------
        evi_data : pandas.DataFrame
            DataFrame with LAI values converted from EVI and measurements
        interpolation_method : str
            'linear' or 'spline' interpolation method for filename suffix
        """
        print(f"\nCreating LAI comparison plot using {interpolation_method} interpolation...")
        
        plots_dir = os.path.join(self.output_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        
        # Set up the plot with white background
        plt.style.use('default')
        fig, ax = plt.subplots(1, 1, figsize=(16, 10))
        
        # Set white background
        ax.set_facecolor('white')
        fig.patch.set_facecolor('white')
        
        dates = evi_data['date']
        
        # Color scheme
        sentinel_color = '#E69F00'  # Orange
        uav_color = '#56B4E9'       # Blue
        field_color = '#000000'     # Black
        lab_color = '#7E2F8E'       # Purple
        
        # Plot Sentinel-2 LAI with 95% prediction intervals
        if 'sentinel_lai_mean' in evi_data.columns and 'sentinel_lai_std' in evi_data.columns:
            sentinel_means = evi_data['sentinel_lai_mean']
            sentinel_stds = evi_data['sentinel_lai_std']
            
            # Filter out NaN values
            valid_mask = ~(np.isnan(sentinel_means) | np.isnan(sentinel_stds))
            if valid_mask.any():
                valid_dates = dates[valid_mask]
                valid_means = sentinel_means[valid_mask]
                valid_stds = sentinel_stds[valid_mask]
                
                print(f"  Plotting Sentinel-2 LAI: {len(valid_dates)} points")
                print(f"    Date range: {valid_dates.min()} to {valid_dates.max()}")
                print(f"    LAI range: {valid_means.min():.3f} to {valid_means.max():.3f}")
                
                # Calculate 95% prediction intervals (±1.96σ for 95% confidence)
                margin_95 = 1.96 * valid_stds
                
                # Plot 95% prediction interval (light colored area) - plot this first so it's behind the line
                '''
                ax.fill_between(valid_dates, 
                               valid_means - margin_95, 
                               valid_means + margin_95, 
                               alpha=0.3, color=sentinel_color, 
                               label='Sentinel-2 LAI (95% PI)')
                '''
                # Plot mean line on top
                ax.plot(valid_dates, valid_means, 
                       color=sentinel_color, linewidth=4,

                       label='Sentinel-2 LAI',
                       markerfacecolor=sentinel_color, 
                       markeredgecolor='white', markeredgewidth=2)
            else:
                print("  No valid Sentinel-2 LAI data to plot")
        
        # Plot UAV LAI with 95% prediction intervals
        if 'uav_lai_mean' in evi_data.columns and 'uav_lai_std' in evi_data.columns:
            uav_means = evi_data['uav_lai_mean']
            uav_stds = evi_data['uav_lai_std']
            
            print(f"  Debug UAV data: {len(uav_means)} total points")
            print(f"    UAV means: {uav_means.describe()}")
            print(f"    UAV stds: {uav_stds.describe()}")
            
            # Filter out NaN values - be more lenient with UAV data
            valid_mask = ~np.isnan(uav_means)  # Only check means, not stds
            if valid_mask.any():
                valid_dates = dates[valid_mask]
                valid_means = uav_means[valid_mask]
                valid_stds = uav_stds[valid_mask]
                
                # Handle NaN stds by setting them to 0
                valid_stds = np.where(np.isnan(valid_stds), 0, valid_stds)
                
                print(f"  Plotting UAV LAI: {len(valid_dates)} points")
                print(f"    Date range: {valid_dates.min()} to {valid_dates.max()}")
                print(f"    LAI range: {valid_means.min():.3f} to {valid_means.max():.3f}")
                
                # Calculate 95% prediction intervals (±1.96σ for 95% confidence)
                margin_95 = 1.96 * valid_stds
                
                # Plot 95% prediction interval (light colored area) - plot this first so it's behind the line
                '''
                ax.fill_between(valid_dates, 
                               valid_means - margin_95, 
                               valid_means + margin_95, 
                               alpha=0.3, color=uav_color, 
                               label='UAV LAI (95% PI)')
                '''
                # Plot mean line on top with more prominent styling
                ax.plot(valid_dates, valid_means, 
                       color=uav_color, linewidth=4, 

                       label='UAV LAI',
                       markerfacecolor=uav_color, 
                       markeredgecolor='white', markeredgewidth=2,
                       zorder=10)  # Ensure it's on top
            else:
                print("  No valid UAV LAI data to plot")
        else:
            print("  UAV LAI columns not found in data")
        

        
        # Plot field measurements
        field_vals = evi_data['field_lai']
        field_mask = ~np.isnan(field_vals)
        if field_mask.any():
            field_dates = dates[field_mask]
            field_values = field_vals[field_mask]
            print(f"  Plotting Field LAI: {len(field_values)} points")
            ax.plot(field_dates, field_values, 
                   color=field_color, linewidth=4,

                   label='Field LAI',
                   markerfacecolor=field_color, 
                   markeredgecolor='white', markeredgewidth=2)
        else:
            print("  No valid Field LAI data to plot")
        
        # Plot lab measurements
        lab_vals = evi_data['lab_lai']
        lab_mask = ~np.isnan(lab_vals)
        if lab_mask.any():
            lab_dates = dates[lab_mask]
            lab_values = lab_vals[lab_mask]
            print(f"  Plotting Lab LAI: {len(lab_values)} points")
            ax.plot(lab_dates, lab_values, 
                   color=lab_color, linewidth=4,

                   label='Lab LAI', 
                   markerfacecolor=lab_color, 
                   markeredgecolor='white', markeredgewidth=2)
        else:
            print("  No valid Lab LAI data to plot")
        
        # Customize the plot
        ax.set_xlabel('Date', fontsize=26)
        ax.set_ylabel('LAI', fontsize=26)
        # Remove title as requested
        
        # Format x-axis dates as DD.MM.YY
        ax.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%d.%m.%y'))
        ax.tick_params(axis='x', rotation=0, labelsize=24)
        ax.tick_params(axis='y', labelsize=24)
        
        # Add legend
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(fontsize=24, frameon=True, fancybox=True, shadow=True,
                  loc="lower center", ncol=4,
                  bbox_to_anchor=(0.5, -0.2))

        # Add grid
        ax.grid(True, alpha=0.3)
        
        # Keep spines (borders) for white background
        for spine in ax.spines.values():
            spine.set_color('black')
            spine.set_linewidth(0)
        
        plt.tight_layout()
        
        # Save the plot
        output_path = os.path.join(plots_dir, f'evi_comparison_plot_{interpolation_method}.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        plt.close()
        
        print(f"✓ Created LAI comparison plot: {output_path}")
        
        return output_path
    
    def create_correlation_analysis(self, evi_data, interpolation_method='spline', correlation_method='pearson'):
        """
        Create correlation analysis between Sentinel LAI and Lab/Field LAI using interpolated values
        
        Parameters:
        -----------
        evi_data : pandas.DataFrame
            DataFrame with interpolated LAI values and measurements
        interpolation_method : str
            'linear' or 'spline' interpolation method for correlation analysis
        correlation_method : str
            'pearson' or 'spearman' correlation method
            
        Returns:
        --------
        dict
            Dictionary with correlation results
        """
        print(f"\nCreating {correlation_method} correlation analysis using {interpolation_method} interpolated values...")
        
        # Prepare data for correlation analysis
        correlation_results = {}
        
        # Use specified interpolation method for each correlation pair
        def compute_corr(col_x, col_y, label):
            if interpolation_method == 'spline':
                x_series = self._spline_interpolate(evi_data['date'], evi_data[col_x])
                y_series = self._spline_interpolate(evi_data['date'], evi_data[col_y])
            else:  # linear
                x_series = self._linear_interpolate(evi_data['date'], evi_data[col_x])
                y_series = self._linear_interpolate(evi_data['date'], evi_data[col_y])
            if x_series is None or y_series is None:
                return None
            # Align indexes
            common_idx = x_series.index.intersection(y_series.index)
            if len(common_idx) < 3:
                return None
            x_vals = x_series.loc[common_idx]
            y_vals = y_series.loc[common_idx]
            # Drop NaNs just in case
            valid = ~(np.isnan(x_vals) | np.isnan(y_vals))
            if valid.sum() < 3:
                return None
            
            # Choose correlation method
            if correlation_method == 'spearman':
                corr, p_val = spearmanr(x_vals[valid], y_vals[valid])
            else:  # pearson (default)
                corr, p_val = pearsonr(x_vals[valid], y_vals[valid])
            
            return {
                'correlation': corr,
                'r_squared': corr**2,
                'p_value': p_val,
                'p_value_fmt': f"{p_val:.3f}" if p_val >= 0.001 else "<0.001",
                'n_points': valid.sum(),
                'significant': p_val < 0.05,
                'method': correlation_method
            }

        pairs = [
            ('sentinel_lai_mean', 'field_lai', 'Sentinel_vs_Field'),
            ('sentinel_lai_mean', 'lab_lai', 'Sentinel_vs_Lab'),
            ('uav_lai_mean', 'field_lai', 'UAV_vs_Field'),
            ('uav_lai_mean', 'lab_lai', 'UAV_vs_Lab'),
            ('sentinel_lai_mean', 'uav_lai_mean', 'Sentinel_vs_UAV'),
            ('field_lai', 'lab_lai', 'Field_vs_Lab'),
        ]
        
        for col_x, col_y, label in pairs:
            if col_x in evi_data.columns and col_y in evi_data.columns:
                res = compute_corr(col_x, col_y, label)
                if res:
                    correlation_results[label] = res
                    print(f"  {label}: corr={res['correlation']:.3f}, n={res['n_points']}")
        
        # Save correlation results
        correlation_df = pd.DataFrame(correlation_results).T
        correlation_output_path = os.path.join(self.output_dir, f'correlation_analysis_{interpolation_method}_{correlation_method}.csv')
        correlation_df.to_csv(correlation_output_path)
        print(f"✓ Saved {correlation_method} correlation analysis: {correlation_output_path}")
        
        # Create correlation scatter plots
        self.create_correlation_plots(evi_data, correlation_results, interpolation_method, correlation_method)
        
        return correlation_results
    
    def create_correlation_plots(self, evi_data, correlation_results, interpolation_method='spline', correlation_method='pearson'):
        """
        Create scatter plots for correlation analysis using interpolated values
        
        Generates scatter plots using daily interpolated LAI values for both X and Y variables
        so the visualisation matches the statistics (n values) used in the correlation computation.
        
        Parameters:
        -----------
        interpolation_method : str
            'linear' or 'spline' interpolation method
        correlation_method : str
            'pearson' or 'spearman' correlation method
        """
        print(f"  Creating {correlation_method} correlation scatter plots with {interpolation_method} interpolated data...")
        
        plots_dir = os.path.join(self.output_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        
        # Mapping of correlations to column pairs and axis labels
        correlation_mappings = {
            'Sentinel_vs_Field':  ('sentinel_lai_mean', 'field_lai',   'Sentinel LAI', 'Field LAI'),
            'Sentinel_vs_Lab':   ('sentinel_lai_mean', 'lab_lai',      'Sentinel LAI', 'Lab LAI'),
            'UAV_vs_Field':      ('uav_lai_mean',      'field_lai',    'UAV LAI',      'Field LAI'),
            'UAV_vs_Lab':        ('uav_lai_mean',      'lab_lai',      'UAV LAI',      'Lab LAI'),
            'Sentinel_vs_UAV':   ('sentinel_lai_mean', 'uav_lai_mean', 'Sentinel LAI', 'UAV LAI'),
            'Field_vs_Lab':      ('field_lai',         'lab_lai',      'Field LAI',    'Lab LAI'),
        }
        
        available_correlations = [c for c in correlation_results.keys() if c in correlation_mappings]
        if not available_correlations:
            print("  No correlation results available to plot – skipping plot generation")
            return
        
        # Determine grid size
        n_plots = len(available_correlations)
        n_cols = 3 if n_plots > 2 else n_plots
        n_rows = int(np.ceil(n_plots / n_cols))
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(7*n_cols, 6*n_rows))
        axes = np.array(axes).flatten()  # ensure 1-d list
        
        # colour palette
        palette = ['#E69F00', '#56B4E9', '#009E73', '#CC79A7', '#F0E442', '#0072B2']
        
        for idx, corr_key in enumerate(available_correlations):
            ax = axes[idx]
            x_col, y_col, x_label, y_label = correlation_mappings[corr_key]
            
            # Interpolate both series to daily frequency using specified method
            if interpolation_method == 'spline':
                x_series = self._spline_interpolate(evi_data['date'], evi_data[x_col])
                y_series = self._spline_interpolate(evi_data['date'], evi_data[y_col])
            else:  # linear
                x_series = self._linear_interpolate(evi_data['date'], evi_data[x_col])
                y_series = self._linear_interpolate(evi_data['date'], evi_data[y_col])
            if x_series is None or y_series is None:
                ax.set_visible(False)
                continue
            common_idx = x_series.index.intersection(y_series.index)
            if len(common_idx) < 3:
                ax.set_visible(False)
                continue
            x_vals = x_series.loc[common_idx]
            y_vals = y_series.loc[common_idx]
            valid = ~(np.isnan(x_vals) | np.isnan(y_vals))
            if valid.sum() < 3:
                ax.set_visible(False)
                continue
            x_vals = x_vals[valid]
            y_vals = y_vals[valid]
            
            ax.scatter(x_vals, y_vals, s=30, alpha=0.6,
                       color=palette[idx % len(palette)], edgecolors='white', linewidth=0.5)
            # 1:1 line
            min_val = min(x_vals.min(), y_vals.min())
            max_val = max(x_vals.max(), y_vals.max())
            ax.plot([min_val, max_val], [min_val, max_val], 'k--', linewidth=1, alpha=0.5)
            
            # Regression line
            if valid.sum() > 3:
                z = np.polyfit(x_vals, y_vals, 1)
                ax.plot([min_val, max_val], np.poly1d(z)([min_val, max_val]),
                        color='red', linewidth=1.5)
                eq_txt = f'y = {z[0]:.3f}x + {z[1]:.3f}'
                ax.text(0.05, 0.92, eq_txt, transform=ax.transAxes, fontsize=9,
                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
            
            # Labels and title
            stats = correlation_results[corr_key]
            if correlation_method == 'spearman':
                ax.set_title(f"{x_label} vs {y_label}\nρ²={stats['r_squared']:.3f}, p={stats['p_value_fmt']} (Spearman)",
                             fontsize=10)
            else:
                ax.set_title(f"{x_label} vs {y_label}\nR²={stats['r_squared']:.3f}, p={stats['p_value_fmt']} (Pearson)",
                             fontsize=10)
            ax.set_xlabel(x_label)
            ax.set_ylabel(y_label)
            ax.grid(alpha=0.3)
            ax.set_facecolor('white')
        
        # Hide unused axes
        for j in range(idx+1, len(axes)):
            axes[j].set_visible(False)
        
        fig.tight_layout()
        plt.savefig(os.path.join(plots_dir, f'correlation_scatter_plots_{interpolation_method}_{correlation_method}.png'),
                    dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        print(f'  ✓ Created correlation_scatter_plots_{interpolation_method}_{correlation_method}.png')
    
    def save_evi_data(self, evi_data, interpolation_method='spline'):
        """
        Save LAI data to CSV
        
        Parameters:
        -----------
        evi_data : pandas.DataFrame
            DataFrame with LAI values converted from EVI and measurements
        interpolation_method : str
            'linear' or 'spline' interpolation method for filename suffix
        """
        # Ensure dates are properly formatted before saving
        evi_data_copy = evi_data.copy()
        evi_data_copy['date'] = pd.to_datetime(evi_data_copy['date']).dt.floor('D')
        
        # Save complete comparison data
        output_path = os.path.join(self.output_dir, 'evi_data', f'evi_comparison_data_{interpolation_method}.csv')
        evi_data_copy.to_csv(output_path, index=False)
        print(f"✓ Saved LAI comparison data: {output_path}")
        
        # Also save final processed LAI data separately
        lai_output_path = os.path.join(self.output_dir, 'lai_data', f'final_lai_comparison_{interpolation_method}.csv')
        evi_data_copy.to_csv(lai_output_path, index=False)
        print(f"✓ Saved final LAI comparison: {lai_output_path}")
    
    def run_complete_analysis(self, plot_only=False, interpolation_method='spline', correlation_method='pearson'):
        """
        Run the complete EVI analysis
        
        Parameters:
        -----------
        plot_only : bool
            If True, only create plots from existing data without reprocessing
        interpolation_method : str
            'linear' or 'spline' interpolation method
        correlation_method : str
            'pearson' or 'spearman' correlation method
        """
        print("EVI Comparison Analysis")
        print("=" * 60)
        print(f"Using {interpolation_method} interpolation method")
        print(f"Using {correlation_method} correlation method")
        print("=" * 60)
        
        if plot_only:
            print("📊 Plot-only mode: Loading existing data and creating plots...")
            return self.run_plot_only_analysis(interpolation_method, correlation_method)
        
        # Step 1: Process Sentinel-2 images
        sentinel_results = self.process_sentinel_images()
        
        # Step 2: Process UAV images
        uav_results = self.process_uav_images()
        
        if sentinel_results is None and uav_results is None:
            print("❌ Error: No data to analyze")
            return None
        
        # Step 3: Interpolate to field dates
        evi_data = self.interpolate_to_field_dates(sentinel_results, uav_results, interpolation_method)
        
        # Step 4: Create comparison plot
        plot_path = self.create_evi_comparison_plot(evi_data, interpolation_method)
        
        # Step 5: Create correlation analysis
        correlation_results = self.create_correlation_analysis(evi_data, interpolation_method, correlation_method)
        
        # Step 6: Save data
        self.save_evi_data(evi_data, interpolation_method)
        
        print(f"\n✅ LAI comparison analysis completed successfully using {interpolation_method} interpolation and {correlation_method} correlation!")
        
        return evi_data, correlation_results
    
    def run_plot_only_analysis(self, interpolation_method='spline', correlation_method='pearson'):
        """
        Run analysis in plot-only mode using existing cached data
        
        Parameters:
        -----------
        interpolation_method : str
            'linear' or 'spline' interpolation method
        correlation_method : str
            'pearson' or 'spearman' correlation method
        
        Returns:
        --------
        tuple
            (evi_data, correlation_results) or None if no cached data available
        """
        print("Loading cached data for plot generation...")
        
        # Try to load cached data
        sentinel_results = self.load_platform_data('sentinel')
        uav_results = self.load_platform_data('uav')
        
        if sentinel_results is None and uav_results is None:
            print("❌ Error: No cached data found. Please run full analysis first.")
            return None
        
        # Check if processed data already exists
        existing_data_path = os.path.join(self.output_dir, 'evi_data', f'evi_comparison_data_{interpolation_method}.csv')
        if os.path.exists(existing_data_path):
            print("📂 Loading existing processed data...")
            try:
                evi_data = pd.read_csv(existing_data_path)
                evi_data['date'] = pd.to_datetime(evi_data['date'])
                
                print(f"✅ Loaded existing data: {len(evi_data)} records")
                print(f"   Date range: {evi_data['date'].min()} to {evi_data['date'].max()}")
                
                # Create plots from existing data
                print("\n📊 Creating plots from existing data...")
                plot_path = self.create_evi_comparison_plot(evi_data, interpolation_method)
                correlation_results = self.create_correlation_analysis(evi_data, interpolation_method, correlation_method)
                
                print("\n✅ Plot-only analysis completed successfully!")
                return evi_data, correlation_results
                
            except Exception as e:
                print(f"⚠️ Error loading existing data: {e}")
                print("Falling back to processing cached raw data...")
        
        # If no processed data exists, process from cached raw data
        print("📊 Processing cached raw data...")
        
        # Step 1: Interpolate to field dates using cached data
        evi_data = self.interpolate_to_field_dates(sentinel_results, uav_results, interpolation_method)
        
        # Step 2: Create comparison plot
        plot_path = self.create_evi_comparison_plot(evi_data, interpolation_method)
        
        # Step 3: Create correlation analysis
        correlation_results = self.create_correlation_analysis(evi_data, interpolation_method, correlation_method)
        
        # Step 4: Save data (for future plot-only runs)
        self.save_evi_data(evi_data, interpolation_method)
        
        print("\n✅ Plot-only analysis completed successfully!")
        
        return evi_data, correlation_results

    def _spline_interpolate(self, dates, values, freq='D'):
        """Return cubic spline interpolated pandas Series at given frequency (default daily)."""
        temp_df = pd.DataFrame({'date': dates, 'val': values}).dropna()
        if len(temp_df) < 3:
            return None  # Need at least 3 points for cubic spline
        t0 = temp_df['date'].min()
        xs = (temp_df['date'] - t0).dt.total_seconds()
        ys = temp_df['val']
        cs = CubicSpline(xs, ys, extrapolate=False)
        new_dates = pd.date_range(t0, temp_df['date'].max(), freq=freq)
        new_xs = (new_dates - t0).total_seconds()
        new_vals = cs(new_xs)
        return pd.Series(new_vals, index=new_dates)
    
    def _linear_interpolate(self, dates, values, freq='D'):
        """Return linear interpolated pandas Series at given frequency (default daily)."""
        temp_df = pd.DataFrame({'date': dates, 'val': values}).dropna()
        if len(temp_df) < 2:
            return None  # Need at least 2 points for linear interpolation
        t0 = temp_df['date'].min()
        xs = (temp_df['date'] - t0).dt.total_seconds()
        ys = temp_df['val']
        interp_func = interp1d(xs, ys, kind='linear', fill_value='extrapolate', bounds_error=False)
        new_dates = pd.date_range(t0, temp_df['date'].max(), freq=freq)
        new_xs = (new_dates - t0).total_seconds()
        new_vals = interp_func(new_xs)
        return pd.Series(new_vals, index=new_dates)

def main(plot_only=False, clear_cache=False):
    """Main function to run the LAI analysis (from EVI conversion) with both interpolation methods"""
    # Configuration
    sentinel_input_dir = '../sentinel2_red_edge_exports'
    uav_input_dir = '../uav'
    shapes_file = '../esu_shapes.geojson'
    output_dir = 'evi_comparison_results'
    
    # Initialize analyzer
    analyzer = EVIAnalyzer(sentinel_input_dir, uav_input_dir, shapes_file, output_dir)
    
    # Clear cache if requested (needed when fixing polygon processing)
    if clear_cache:
        print("🗑️ Clearing cache to ensure polygon fix takes effect...")
        analyzer.clear_cache()
    
    results = {}
    
    # Run analysis with both interpolation methods and both correlation methods
    for interp_method in ['linear', 'spline']:
        for corr_method in ['pearson', 'spearman']:
            method_key = f"{interp_method}_{corr_method}"
            print(f"\n{'='*80}")
            print(f"RUNNING ANALYSIS WITH {interp_method.upper()} INTERPOLATION + {corr_method.upper()} CORRELATION")
            print(f"{'='*80}")
            
            # Run complete analysis
            evi_data, correlation_results = analyzer.run_complete_analysis(
                plot_only=plot_only, 
                interpolation_method=interp_method,
                correlation_method=corr_method
            )
            
            results[method_key] = {
                'evi_data': evi_data,
                'correlation_results': correlation_results
            }
    
    print(f"\n{'='*80}")
    print("ANALYSIS COMPLETE - All variants generated:")
    print("📊 Linear interpolation + Pearson correlation")
    print("📊 Linear interpolation + Spearman correlation")
    print("📊 Spline interpolation + Pearson correlation")
    print("📊 Spline interpolation + Spearman correlation")
    print("Check the output directory for method-specific files")
    print(f"{'='*80}")
    
    return results

if __name__ == "__main__":
    # Clear cache to ensure polygon fix takes effect
    main(plot_only=False)
