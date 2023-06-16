// Data0: Double check the time series of specific coordinates

var pt = ee.Geometry.Point([97.16455, 21.25997]);
Map.centerObject(pt, 11);

var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterDate('2013-01-01', '2020-02-01')
    .filterBounds(pt);

//l8 = l8.map(applyScaleFactors);
//l8 = l8.filter(ee.Filter.lte('CLOUD_COVER_LAND', 1))

function mask_landsat_sr(image) {
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
  var cloudShadowBitMask = (1 << 4);
  var cloudsBitMask = (1 << 3);
  // Get the pixel QA band.
  var qa = image.select('QA_PIXEL');
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  var saturationMask = image.select('QA_RADSAT').eq(0);
  return image.updateMask(mask).updateMask(saturationMask);
}

var l8 = l8.map(mask_landsat_sr).select('SR.*')

var l8NDVI = l8.map(function(image) {
  // clip image
  //var clipped = image.clip(pt)
  
  var nir = image.select(['SR_B5']).multiply(0.0000275).add(-0.2);
  var red = image.select(['SR_B4']).multiply(0.0000275).add(-0.2);
  
  // Compute NDVI.
  var ndvi = image.expression(
                '(nir - red) / (nir + red)', {
                  'nir': nir,
                  'red': red
                });

  // Return the masked image with an NDVI band.
  var out = image.addBands(ndvi);
  var outTrue = out.select("SR_B5_1").rename("NDVI")
  return outTrue
});

//////////////////////////////////////////////////////////////
// Sentinel 2
//////////////////////////////////////////////////////////////
// Source: https://developers.google.com/earth-engine/tutorials/community/sentinel-2-s2cloudless

var AOI = pt
var START_DATE = '2015-06-23'
var END_DATE = '2020-02-01'
var CLOUD_FILTER = 50 //Maximum image cloud cover percent allowed in image collection
var CLD_PRB_THRESH = 30 //Cloud probability (%); values greater than are considered cloud
var NIR_DRK_THRESH = 0.15 //Near-infrared reflectance; values less than are considered potential cloud shadow
var CLD_PRJ_DIST = 2 //Maximum distance (km) to search for cloud shadows from cloud edges
var BUFFER = 50 //Distance (m) to dilate the edge of cloud-identified objects

var get_s2_sr_cld_col = function(aoi, start_date, end_date){
    // Import and filter S2 SR.
    var s2_sr_col = (ee.ImageCollection('COPERNICUS/S2') //S2_SR
        .filterBounds(aoi)
        .filterDate(start_date, end_date))
        //.filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', CLOUD_FILTER)))
    //print(s2_sr_col)

    // Import and filter s2cloudless.
    var s2_cloudless_col = (ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
        .filterBounds(aoi)
        .filterDate(start_date, end_date))
    //print(s2_cloudless_col)

    // Join the filtered s2cloudless collection to the SR collection by the 'system:index' property.
    return ee.ImageCollection(ee.Join.saveFirst('s2cloudless').apply({
        'primary': s2_sr_col,
        'secondary': s2_cloudless_col,
        'condition': ee.Filter.equals({
            'leftField': 'system:index',
            'rightField': 'system:index'
        })
    }))
}
var add_cloud_bands = function(img){
    // Get s2cloudless image, subset the probability band.
    var cld_prb = ee.Image(img.get('s2cloudless')).select('probability')

    // Condition s2cloudless by the probability threshold value.
    var is_cloud = cld_prb.gt(CLD_PRB_THRESH).rename('clouds')

    // Add the cloud probability layer and cloud mask as image bands.
    return img.addBands(ee.Image([cld_prb, is_cloud]))
}
var add_shadow_bands = function(img){
    // Identify water pixels from the SCL band.
    //var not_water = img.select('SCL').neq(6)

    // Identify dark NIR pixels that are not water (potential cloud shadow pixels).
    var SR_BAND_SCALE = 1e4
    //var dark_pixels = img.select('B8').lt(NIR_DRK_THRESH*SR_BAND_SCALE).multiply(not_water).rename('dark_pixels')
    var dark_pixels = img.select('B8').lt(NIR_DRK_THRESH*SR_BAND_SCALE).rename('dark_pixels')

    // Determine the direction to project cloud shadow from clouds (assumes UTM projection).
    var shadow_azimuth = ee.Number(90).subtract(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')));

    // Project shadows from clouds for the distance specified by the CLD_PRJ_DIST input.
    var cld_proj = (img.select('clouds').directionalDistanceTransform(shadow_azimuth, CLD_PRJ_DIST*10)
        .reproject({'crs': img.select(0).projection(), 'scale': 100})
        .select('distance')
        .mask()
        .rename('cloud_transform'))

    // Identify the intersection of dark pixels with cloud shadow projection.
    var shadows = cld_proj.multiply(dark_pixels).rename('shadows')

    // Add dark pixels, cloud projection, and identified shadows as image bands.
    return img.addBands(ee.Image([dark_pixels, cld_proj, shadows]))
}
var add_cld_shdw_mask = function(img){
    // Add cloud component bands.
    var img_cloud = add_cloud_bands(img)

    // Add cloud shadow component bands.
    var img_cloud_shadow = add_shadow_bands(img_cloud)

    // Combine cloud and shadow mask, set cloud and shadow as value 1, else 0.
    var is_cld_shdw = img_cloud_shadow.select('clouds').add(img_cloud_shadow.select('shadows')).gt(0)
    //var is_cld = img_cloud.select('clouds').gt(0).rename('cloudmask')

    // Remove small cloud-shadow patches and dilate remaining pixels by BUFFER input.
    // 20 m scale is for speed, and assumes clouds don't require 10 m precision.
    
    var is_cld_shdw = (is_cld_shdw.focal_min(2).focal_max(BUFFER*2/20)
        .reproject({'crs': img.select([0]).projection(), 'scale': 20})
        .rename('cloudmask'))
    
    
    // Add the final cloud-shadow mask to the image.
    return img_cloud_shadow.addBands(is_cld_shdw)
    //return img_cloud.addBands(is_cld)
}
var apply_cld_shdw_mask = function(img){
    // Subset the cloudmask band and invert it so clouds/shadow are 0, else 1.
    var not_cld_shdw = img.select('cloudmask').not()

    // Subset reflectance bands and update their masks, return the result.
    return img.select('B.*').updateMask(not_cld_shdw)
}


var s2_sr_cld_col_eval = get_s2_sr_cld_col(AOI, START_DATE, END_DATE)
var s2Mask = s2_sr_cld_col_eval.map(add_cld_shdw_mask)
                             .map(apply_cld_shdw_mask)


/*
var tstimg = s2_sr_cld_col_eval.first()
print(tstimg)
print(ee.Image(tstimg.get('s2cloudless')).select('probability'))
*/


// Function to calculate and add an NDVI band
var addNDVI = function(image) {
return image.addBands(image.normalizedDifference(['B8A', 'B4']));
};

// Add NDVI band to image collection
var s2ndvi = s2Mask.map(addNDVI);

print(ui.Chart.image.series({
  imageCollection: l8NDVI.select('NDVI'),
  region: pt,
  reducer: ee.Reducer.first(),
  scale: 30
}).setOptions({title: 'L8: Cloud-masked NDVI'}));

print(ui.Chart.image.series({
  imageCollection: s2ndvi.select(['nd']),
  region: pt,
  reducer: ee.Reducer.first(),
  scale: 10
}).setOptions({title: 'S2: Cloud-masked NDVI'}));

var visualization = {
  bands: ['SR_B4', 'SR_B3', 'SR_B2'],
  min: 0.0,
  max: 0.3,
};
