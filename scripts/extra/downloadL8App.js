//GEE: use 1 extent to check. 
//Get L8 C2, crop, export

//NOTE that the input coords here are already buffered by 1000m from the actual extent
// we are using for the landscape analysis
var xmin = 97.1388448997317
var xmax = 97.4192548921057
var ymin = 24.7835951329582
var ymax = 24.9837340984572

var polygon = ee.Geometry.Polygon([
  [xmin, ymin], [xmin, ymax],
  [xmax, ymax], [xmax, ymin]
]);

Map.addLayer(polygon);
Map.centerObject(polygon, 11);

var col2019 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
                .filterDate("2018-12-15", "2019-01-15")
                .filterBounds(polygon)

// A function that scales and masks Landsat 8 (C2) surface reflectance images.
function prepSrL8(image) {
  // Develop masks for unwanted pixels (fill, cloud, cloud shadow).
  var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111', 2)).eq(0);
  var saturationMask = image.select('QA_RADSAT').eq(0);

  // Apply the scaling factors to the appropriate bands.
  var getFactorImg = function(factorNames) {
    var factorList = image.toDictionary().select(factorNames).values();
    return ee.Image.constant(factorList);
  };
  var scaleImg = getFactorImg([
    'REFLECTANCE_MULT_BAND_.|TEMPERATURE_MULT_BAND_ST_B10']);
  var offsetImg = getFactorImg([
    'REFLECTANCE_ADD_BAND_.|TEMPERATURE_ADD_BAND_ST_B10']);
  var scaled = image.select('SR_B.|ST_B10').multiply(scaleImg).add(offsetImg);

  // Replace original bands with scaled bands and apply masks.
  return image.addBands(scaled, null, true)
    .updateMask(qaMask).updateMask(saturationMask);
}

// Landsat 8 Collection 2 surface reflectance images of interest.
var colBefore = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(polygon)
  .filterDate('2018-07-01', '2019-01-01')
  .map(prepSrL8)
  .select('SR.*')
  .median();
  
var colAfter = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(polygon)
  .filterDate('2019-05-01', '2020-07-01')
  .map(prepSrL8)
  .select('SR.*')
  .median();

Map.addLayer(colBefore.clip(polygon), {bands: ['SR_B4', 'SR_B3', 'SR_B2'], min: 0, max: 0.1});
Map.addLayer(colAfter.clip(polygon), {bands: ['SR_B4', 'SR_B3', 'SR_B2'], min: 0, max: 0.1});

// Actual main downloading portion of the script
var collection = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
                .filterDate("2015-06-23", "2020-02-01")
                .filterBounds(polygon)
                .select("SR_B4", "SR_B5", "QA_PIXEL", "QA_RADSAT", "SR_QA_AEROSOL")

var clipped = collection.map(function(img){
  return img.clip(polygon)
})
//print(clipped)

var visParams = {
  bands: ['SR_B4', 'SR_B5', 'QA_PIXEL'],
  min: 0,
  max: 0.1
};
//Map.addLayer(clipped.first(), visParams);


//var batch = require('users/fitoprincipe/geetools:batch');
//batch.ImageCollection.toDrive(clipped, 'imagesL8app', {scale: 30, type: 'uint16'})
  
var tools = require('users/fitoprincipe/geetools:tools');    
/*var downloadAll = function(collection, folder, options) {
    
    var defaults = {
      scale: 1000,
      maxPixels: 1e13,
      type: 'float',
      region: null,
      name: null,
    }
    
    var opt = tools.get_options(defaults, options)
    
    var n = collection.size().getInfo();
    
    var colList = collection.toList(n);
    
    var colID = opt.name || collection.getInfo()['id'] || ""
    colID = colID.replace('/','_')
    
    for (var i = 0; i < n; i++) {
      var img = ee.Image(colList.get(i));
      var id = img.id().getInfo() || colID+'_image_'+i.toString();
      var region = opt.region || img.geometry().bounds().getInfo()["coordinates"];
      
      //var imtype = IMAGE_TYPES(img, opt.type)
      var imgtype = img.uint16()
      
      Export.image.toDrive({
        image: imgtype,
        description: id,
        folder: folder,
        fileNamePrefix: id,
        scale: opt.scale,
        crs: "EPSG: 4326",
        maxPixels: opt.maxPixels})
    }
  }*/

/*
Export.image.toDrive({
        image: clipped.first(),
        folder: 'imagesL8app',
        scale: 30,
        crs: 'EPSG: 4326',
        maxPixels: 1e13})
*/
downloadAll(clipped, 'imagesL8app', {scale: 30})