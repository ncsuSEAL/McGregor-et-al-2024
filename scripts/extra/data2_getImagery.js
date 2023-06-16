// PURPOSE: extract time series of vegetation indices for training locations, as a preliminary
// model to inform further dissertation studies
// SOURCE: Ian McGregor, NCSU
// FIRST WRITTEN: March-April 2020, last updated August 2021

//BEFORE RUNNING//
// Upload trainingDataPoints.csv as an asset to GEE, and then add to this script as
// the variable "table"

//METHODS//
//we have a series of training points.
//we want to get a full NDVI timeseries for each of those points for with a 1km buffer
//the time series should only be for the lifespan of sentinel 2
print('original csv', table);

//---------------------------------------------------------------------------//
// Step 0: Define the function (see below at bottom)
//---------------------------------------------------------------------------//
var data = function(satellite, filename, fullbandnames, pointGroup, centroid, country){
  
  ///---------------------------------------------------------------------------///
  /// Step 1: Define the ROI for map display and shapefiles
  ///---------------------------------------------------------------------------///
  //Note I was using this earlier but now I'm just doing centerObject over the points
  var display = ee.Geometry.Point(centroid).buffer(300000)
  
  print("Focal country: ", country)
  
  var dataset = ee.FeatureCollection('USDOS/LSIB/2013') //country boundaries
    .filterBounds(display)
    .first();

  var image1 = [ee.Image("users/imcgreg/glad19").clip(dataset),
                ee.Feature(dataset)]
  //Map.addLayer(image1, vizParams);
  var vizParams = {opacity: 0.3};
  ///---------------------------------------------------------------------------///
  /// Step 2: Get training points (from asset) and create index from them
  ///---------------------------------------------------------------------------///
  var getMetadata = function(feature) {
      var lon = feature.get('coordX');
      var lat = feature.get('coordY');
      return ee.Feature(ee.Geometry.Point([lon, lat]), {
        'featureID': ee.String(feature.get('pointid'))
      });
    };
  
  //NOTE the first column of the table (before importing) MUST be a random column that is not
  //coordX, coordY, or pointid, otherwise GEE will not recognize it (no idea why)
  
  var points = ee.FeatureCollection(
      ee.FeatureCollection(table)
          .filter(ee.Filter.stringContains('location', country))
    ).map(getMetadata);
  
  /*
  var points = ee.Feature(ee.Geometry.Point([-72.171937, 42.537326]), {
        'featureID': ee.String('harvard')
      })*/
  
  print('points', points.get(0).coordinates); 
  //print(points.first())
  //print(points.first().get('pointid'));
  Map.addLayer(points, {color: "blue"}); 
  
  var chatthin = ee.FeatureCollection('users/imcgreg/chatthin');
  chatthin = chatthin.geometry();
  Map.addLayer(chatthin, {color: 'red'}, 'Chatthin');
  
  var mahamyaing = ee.FeatureCollection('users/imcgreg/mahamyaing');
  mahamyaing = mahamyaing.geometry();
  Map.addLayer(mahamyaing, {color: 'red'}, 'Mahamyaing');

  //id and coordinates for each of the points
  var ids = points.aggregate_array('featureID');
  //var ids = ee.List([0])
  var coords = points.geometry().coordinates();

  print('ids', ids);
  print('coords', coords);
  
  Map.centerObject(points)

  // make index for the points
  //var index = ee.List([0])
  var index = ee.List.sequence(0,coords.size().subtract(1)); print(index);
  print('loop size', index.size());
  
  ///---------------------------------------------------------------------------///
  /// Step 3: Obtain satimgs filtered to all points for L8 and S2 (+S2cloudless)
  ///---------------------------------------------------------------------------///
  // "filename" is defined at the bottom as part of the overall function
  var s1 = filename == 'sentinel1'; print('is s1?', s1);
  var l8 = filename == 'landsat8sr'; print('is l8?', l8);
  
  var date_start = ee.String(
    ee.Algorithms.If(ee.String(country).equals("myanmar"), 
      ee.Algorithms.If(l8, '2013-02-01', '2015-06-23'),
      '2019-01-01'));
  var date_end = ee.String(
    ee.Algorithms.If(ee.String(country).equals("myanmar"), '2020-02-01',
      ee.Algorithms.If(ee.String(country).equals("peru"), '2022-02-01', '2022-04-22')));

  print("Start date: ", date_start); print("End date: ", date_end);

  var satimgs_init = ee.ImageCollection(
    ee.Algorithms.If(s1, null, ee.ImageCollection(satellite)
                  .filterDate(date_start, date_end)
                  .filterBounds(points)));
  
  ///---------------------------------------------------------------------------///
  /// Step 3a: Define helper functions
  ///---------------------------------------------------------------------------///
  // Define the join filter for combining VV and VH / image collections
  var filter = ee.Filter.equals({
    leftField: 'system:time_start',
    rightField: 'system:time_start'
  });

  var mergeByFeature = function(feature) {
    return ee.Image.cat(ee.Image(feature.get('primary')), 
                        ee.Image(feature.get('secondary')));
  };
  
  function renameBands(oldname, newname) {
      var nest = function(image){
        return image.select(oldname).rename(newname);
      };
      return nest;
    }
    
  // properties to carry through the processing
  var save_props = ['system: id', 'system:time_start', 
            'instrumentMode', 'resolution_meters', 
            'resolution', 'orbitProperties_pass']
  
  ///---------------------------------------------------------------------------///
  /// Step 3b: Apply functions to S1 imagery (if satellite = S1)
  ///---------------------------------------------------------------------------///
  //We will apply best practice S1 ARD methods as defined in Mullissa et al 2021
  /*
  // filter to keep VV and VH
  var processS1 = function(s1imgs){
    var filtered = s1imgs.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
              .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
              .filter(ee.Filter.eq('instrumentMode', 'IW'));
    print('filtered', filtered.first());
    
    ///---------------------------------------------------------------------------///
    /// Step 3bi: Apply border noise correction
    ///---------------------------------------------------------------------------///
    /// we do not need to filter by asc/desc because this is accounted for by 
    /// orbit number during multitemporal speckle filtering
    var s1_B = filtered.map(border_noise_correction.f_mask_edges)
    print('S1 border noise', s1_B.first())
     
    ///---------------------------------------------------------------------------///
    /// Step 3bii: Apply speckle filtering
    ///---------------------------------------------------------------------------///
    // This is multitemporal boxcar (Quegan and Yu 2001) + Refined Lee
    var speckleKernel = 7
    var speckleFilter = 'REFINED LEE'
    var speckleImgNum = 10
    
    var s1_BS = ee.ImageCollection(speckle_filter
                .MultiTemporal_Filter(s1_B, speckleKernel, speckleFilter, speckleImgNum))
    
    print('S1 RLF', s1_BS.first())
     
    ///---------------------------------------------------------------------------///
    /// Step 3biii: Apply radiometric terrain normalization (slope correction)
    ///---------------------------------------------------------------------------///
    // From Vollrath et al 2020
    var terrainModel = 'VOLUME'
    var DEM = ee.Image('USGS/SRTMGL1_003')
    var terrainBuffer = 0
    var s1_BSR = ee.ImageCollection(terrain_flattening
                  .slope_correction(s1_BS, terrainModel, DEM, terrainBuffer))
  
    print(s1_BSR.first())
    return(s1_BSR)
  }
  */
  
  ///---------------------------------------------------------------------------///
  /// Step 3b: Define functions for cloudmask for s2 (using s2 cloudless)
  ///---------------------------------------------------------------------------///
  // "filename" is defined at the bottom as part of the overall function
  var s2 = filename == "sentinel2toa" || filename == "sentinel2sr"; 
  print('is s2?', s2);

  var CLD_PRB_THRESH = 30; //Cloud probability (%); values greater than are considered cloud
  var NIR_DRK_THRESH = 0.15; //Near-infrared reflectance; values less than are considered potential cloud shadow
  var CLD_PRJ_DIST = 2; //Maximum distance (km) to search for cloud shadows from cloud edges
  var BUFFER = 50; //Distance (m) to dilate the edge of cloud-identified objects
  
  var get_s2_sr_cld_col = function(s2coll, start, end){
    // Import and filter s2cloudless.
    var s2_cloudless_col = (ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
        .filterBounds(points)
        .filterDate(start, end));
    //print(s2_cloudless_col)

    // Join the filtered s2cloudless collection to the SR collection by the 'system:index' property.
    return ee.ImageCollection(ee.Join.saveFirst('s2cloudless').apply({
        'primary': s2coll,
        'secondary': s2_cloudless_col,
        'condition': ee.Filter.equals({
            'leftField': 'system:index',
            'rightField': 'system:index'
        })
    }));
  };
  
  //add in the cloud bands
  var add_cloud_bands = function(img){
    // Get s2cloudless image, subset the probability band.
    var cld_prb = ee.Image(img.get('s2cloudless')).select('probability');

    // Condition s2cloudless by the probability threshold value.
    var is_cloud = cld_prb.gt(CLD_PRB_THRESH).rename('clouds');

    // Add the cloud probability layer and cloud mask as image bands.
    return img.addBands(ee.Image([cld_prb, is_cloud]));
  };

  // add shadow bands
  var add_shadow_bands = function(img){
    // Identify water pixels from the SCL band.
    //var not_water = img.select('SCL').neq(6) //only if doing SR and have water/land

    // Identify dark NIR pixels that are not water (potential cloud shadow pixels).
    var SR_BAND_SCALE = 1e4;
    //var dark_pixels = img.select('B8').lt(NIR_DRK_THRESH*SR_BAND_SCALE).multiply(not_water).rename('dark_pixels')
    var dark_pixels = img.select('B8').lt(NIR_DRK_THRESH*SR_BAND_SCALE).rename('dark_pixels');

    // Determine the direction to project cloud shadow from clouds (assumes UTM projection).
    var shadow_azimuth = ee.Number(90).subtract(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')));

    // Project shadows from clouds for the distance specified by the CLD_PRJ_DIST input.
    var cld_proj = (img.select('clouds').directionalDistanceTransform(shadow_azimuth, CLD_PRJ_DIST*10)
        .reproject({'crs': img.select(0).projection(), 'scale': 100})
        .select('distance')
        .mask()
        .rename('cloud_transform'));

    // Identify the intersection of dark pixels with cloud shadow projection.
    var shadows = cld_proj.multiply(dark_pixels).rename('shadows');

    // Add dark pixels, cloud projection, and identified shadows as image bands.
    return img.addBands(ee.Image([dark_pixels, cld_proj, shadows]));
  };

  // generate the full cloud mask
  var add_cld_shdw_mask = function(img){
    // Add cloud component bands.
    var img_cloud = add_cloud_bands(img);

    // Add cloud shadow component bands.
    var img_cloud_shadow = add_shadow_bands(img_cloud);

    // Combine cloud and shadow mask, set cloud and shadow as value 1, else 0.
    var is_cld_shdw = img_cloud_shadow.select('clouds').add(img_cloud_shadow.select('shadows')).gt(0);
    //var is_cld = img_cloud.select('clouds').gt(0).rename('cloudmask')

    // Remove small cloud-shadow patches and dilate remaining pixels by BUFFER input.
    // 20 m scale is for speed, and assumes clouds don't require 10 m precision.
    
    var is_cld_shdw = (is_cld_shdw.focal_min(2).focal_max(BUFFER*2/20)
        .reproject({'crs': img.select([0]).projection(), 'scale': 20})
        .rename('cloudmask'));
    
    
    // Add the final cloud-shadow mask to the image.
    //return img_cloud_shadow.addBands(is_cld_shdw) //gives all bands back
    return img_cloud.addBands(is_cld_shdw); //gives only shadow band
  };

  //---------------------------------------------------------------------------///
  // Step 3c: Do S2 preprocessing, otherwise move on
  //---------------------------------------------------------------------------///
  var applyS2Cloud = function(s2coll){
    var s2wCloud = get_s2_sr_cld_col(s2coll, date_start, date_end);
    var s2wCloudMask = s2wCloud.map(add_cld_shdw_mask);
    return s2wCloudMask;
  };
  
  var satimgs = ee.ImageCollection(
    ee.Algorithms.If(s2, applyS2Cloud(satimgs_init), satimgs_init));

  if(s1){print('not S2 or L8')} else {print('testsats', satimgs.first())}
  if(s1){print('not S2 or L8')} else {print('bands', satimgs.first().bandNames())}
  
   ///---------------------------------------------------------------------------///
  /// Step 4: Define variables and functions for S1 processing
  ///---------------------------------------------------------------------------///
  //import modules from Mullissa et al 2021
  var wrapper = require('users/imcgreg/dissertation:s1/s1ARD_Mullissa/s1_wrapper');
  var helper = require('users/adugnagirma/gee_s1_ard:utilities');
  var speckle_filter = require('users/adugnagirma/gee_s1_ard:speckle_filter');
  var terrain_flattening = require('users/adugnagirma/gee_s1_ard:terrain_flattening');
  var border_noise_correction = require('users/adugnagirma/gee_s1_ard:border_noise_correction');
  
  var processS1 = function(geomPoint){
    var parameters = {//1. Data Selection
              START_DATE: date_start,
              STOP_DATE: date_end,
              POLARIZATION:'VVVH',
              ORBIT : 'BOTH',
              GEOMETRY: geomPoint,
              //2. Additional Border noise correction
              APPLY_ADDITIONAL_BORDER_NOISE_CORRECTION: true,
              //3.Speckle filter
              APPLY_SPECKLE_FILTERING: true,
              SPECKLE_FILTER_FRAMEWORK: 'MULTI',
              SPECKLE_FILTER: 'REFINED LEE',
              SPECKLE_FILTER_KERNEL_SIZE: 7,
              SPECKLE_FILTER_NR_OF_IMAGES: 10,
              //4. Radiometric terrain normalization
              APPLY_TERRAIN_FLATTENING: true,
              DEM: ee.Image('USGS/SRTMGL1_003'),
              TERRAIN_FLATTENING_MODEL: 'VOLUME',
              TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER: 0,
              //5. Output
              FORMAT : 'DB', //if LINEAR, then the values are bounded to 0
              CLIP_TO_ROI: false
    };
  //Preprocess the S1 collection
  var s1_preprocess = wrapper.s1_preproc(parameters);
  var s1_preprocess = s1_preprocess[1]; //2 outputs, take the 2nd one

  //Convert format for export
  if (parameters.FORMAT=='DB'){
    s1_preprocess = s1_preprocess.map(helper.lin_to_db);
    }

  return(s1_preprocess)
  //print(s1_preprocess.first())
  }

  ///---------------------------------------------------------------------------///
  /// Step 5: Get band values and properties per point per sensor
  ///---------------------------------------------------------------------------///
  
  // define dates function
  var get_dates = function(collection) {
    return ee.List(collection.toList(collection.size()).map(function(img){
      return ee.Image(img).date().format('YYYY-MM-dd');
    }));
  };
  
  // time conversion functions
  var microSecToHour = function(aNum){
      return aNum.multiply(ee.Number(2.77778).multiply(ee.Number(10).pow(-10)));
  };
  
  var microSecToDays = function(aNum){
    return aNum.multiply(ee.Number(1.15741).multiply(ee.Number(10).pow(-11)));
  };
  
  ////--------------------------------------------------////
  //// Step 5a: Create index to run the loop
  ////--------------------------------------------------////
  
  // first we need to make a list containing the imageColls for each point
  var test1 = ee.List(ee.Algorithms.If(s1, null, 
  index.map(function(q){
    var sub = satimgs.filterBounds(ee.Geometry.Point(coords.get(q)));
    return sub;
  })));
  
  /*
  //run test for first one to inspect elements if need be
  var firstcoll = ee.ImageCollection(test1.get(0));
  var img = firstcoll.first();
  var bands0 = img.bandNames();
  print(img.select(bands0))
  
  print('img', img);
  print('bands', bands0); 
  
  var id11 = ee.String(ids.get(0)); print('id11', id11);
  var date1 = img.date().format('yyyy-MM-dd'); print('date1', date1);
  var satellite1 = ee.String(img.get("system:index")); print('satellite1', satellite1);
  var thispoint1 = ee.Geometry.Point(coords.get(0)); print('thispoint1', thispoint1);
  var spec1 = thispoint1.coordinates(); print('spec1', spec1);
    
  var value = img.select(bands0)
    .reduceRegion(ee.Reducer.first(), thispoint1); print('value', value);
  */
  
  ////--------------------------------------------------////
  //// Step 5b: Run the loop
  ////--------------------------------------------------////
  /*
  var testLoop = index.map(function(ind){
    var thispoint = ee.Geometry.Point(coords.get(ind))
    var coll = ee.ImageCollection(
      ee.Algorithms.If(s1, processS1(thispoint), test1.get(ind))
    );
    return(coll)
  })
  print(ee.ImageCollection(testLoop.get(0)).first())
  kjghjknjhhk
  */
  
  var loop = index.map(function(loopN){
    //first, get the imageColl for each point from either satimgs (S2, L8) or
    ///by applying the ARD functions (S1)
    var thispoint = ee.Geometry.Point(coords.get(loopN))
    var s1coll = ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT')
      .filter(ee.Filter.eq('instrumentMode', 'IW'))
      .filter(ee.Filter.eq('resolution_meters', 10))
      .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
      .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
      .filterDate(date_start, date_end)
      .filterBounds(thispoint);
    
    var coll = ee.ImageCollection(
      ee.Algorithms.If(s1, processS1(thispoint), test1.get(loopN))
    );
    /*
    var coll = ee.ImageCollection(
      ee.Algorithms.If(s1, s1coll, test1.get(loopN))
    );
    */
  
    ////--------------------------------------------------////
    //// Step 5c-i: Retrieve values for each point
    ////--------------------------------------------------////
    
    //get the value for each band for every image in each imagecollection
    // SOME DATES HAVE NULL VALUES thus show up as not having the information
    var getvals = function(img1){
      var bands = img1.bandNames();
      var id1 = ee.String(ids.get(loopN));
      var date = img1.date().format('yyyy-MM-dd');
      var satellite = ee.String(ee.Algorithms.If(s1, 
                                ee.String(img1.get('system:index')),
                                ee.String(img1.get("system:id"))));
      var spec = thispoint.coordinates();
      var value = img1.select(bands)
        .reduceRegion(ee.Reducer.first(), thispoint);
      
      ////--------------------------------------------------////
      //// Step 5c-ii: Bring in S1 properties
      ////--------------------------------------------------////
      var mode = ee.Algorithms.If(s1,
                  ee.String(img1.get('instrumentMode')))
      var res_m = ee.Algorithms.If(s1,
                  ee.String(img1.get('resolution_meters')))
      var res_type = ee.Algorithms.If(s1,
                  ee.String(img1.get('resolution')))
      var dir = ee.Algorithms.If(s1,
                  ee.Number(img1.get('orbitProperties_pass')))
      
      ////--------------------------------------------------////
      //// Step 5c-iii: Convert to FeatureCollection, keep bands
      ////--------------------------------------------------////
      var output = ee.Algorithms.If(s1,
        ee.Feature(ee.Geometry.Point(coords.get(loopN)), 
      {date: date, satellite: satellite, lat: spec.get(0), lon: spec.get(1),
        pointid: id1, mode: mode, res_m: res_m, res_type: res_type, 
        direction: dir})
        .set(value),
      ee.Feature(ee.Geometry.Point(coords.get(loopN)), 
      {date: date, satellite: satellite, lat: spec.get(0), lon: spec.get(1),
        pointid: id1})
        .set(value));
        
    return output
    };
    var newft = coll.map(getvals);
    return(newft);
  });
  print('loopresult', loop)
  //print('loop ex', ee.FeatureCollection(loop.get(2)));
  
  ////------------------------------------------------------------------////
  //// Step 5d: Run the same loop but for observation / ingestion dates
  ////------------------------------------------------------------------////
  
  // see other script ('training_dates')

  ////--------------------------------------------------////
  //// Step 5e: Cast loop to FC, combine (flatten) into one (this is what we export)
  ////--------------------------------------------------////
  var test = ee.FeatureCollection(loop).flatten();

  ///---------------------------------------------------------------------------///
  /// Step 6: Export to csv
  ///---------------------------------------------------------------------------///
  Export.table.toDrive(
    {collection: test,
    description: filename,
    folder: 'data',
    fileNamePrefix: filename,
    selectors: fullbandnames
    }
  );
};

//---------------------------------------------------------------------------//
// Step 0a: Define the list inputs
//---------------------------------------------------------------------------//
var disp = 
  [[95.219623, 23.029892],
   [-70.49, -12.35],
   [-78.728916, 35.824405]];

var countryName = ["myanmar", "peru", "us"];
   
var sat =
  ["COPERNICUS/S1_GRD_FLOAT",   //sentinel 1 SAR-GRD cband (2014-10-03) 4825 features
    "COPERNICUS/S2",            //sentinel 2 TOA Level 1C (2015-06-23) 7425 features
    "COPERNICUS/S2_SR",         //sentinel2 SR Level 2A (2017-03-28) 2810 features
    "LANDSAT/LC08/C01/T1_TOA",  //landsat8 TOA (April 2013) 2324 features
    "LANDSAT/LC08/C02/T1_L2",   //L8 Collection 2 SR
    "LANDSAT/LC09/C02/T1_L2"    //L9 Collection 2 SR
  ];
  
var filenames = 
  ["sentinel1", 
  "sentinel2toa", 
  "sentinel2sr", 
  "landsat8toa", 
  "landsat8sr",
  "landsat9sr"];

var basebands = ['satellite', 'date', 'lat', 'lon', 'pointid'] + ','; //need this last comma!

var bandnm_s1 = basebands + ['VH', 'VV', 'angle', 'mode', 'res_m', 'res_type', 'direction'];
var bandnm_s2toa = basebands + ['B1','B2','B3','B4','B5','B6','B7','B8','B8A',
                                'B9','B10','B11','B12',
                                'QA60', 'cloudmask']; ///QA10/20 are always empty, QA60 is cloud mask
var bandnm_s2sr = basebands + ['B1','B2','B3','B4','B5','B6','B7','B8','B8A',
                                'B9','B11','B12','AOT','WVP','MSK_CLDPRB',
                                'MSK_SNWPRB', 'SCL', 'QA60', 'cloudmask']; //QA10/20 are always empty, QA60 is cloud mask
var bandnm_l8toa = basebands + ['B1','B2','B3','B4','B5','B6','B7','B8','B9',
                                'B10','B11','BQA'];
var bandnm_l8sr = basebands + ['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7','ST_B10',
                              'QA_PIXEL', 'QA_RADSAT', 'SR_QA_AEROSOL'];
var bandnm_l9sr = basebands + ['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7','ST_B10',
                              'QA_PIXEL', 'QA_RADSAT', 'SR_QA_AEROSOL'];

var allbands = [bandnm_s1, 
  bandnm_s2toa, 
  bandnm_s2sr, 
  bandnm_l8toa,
  bandnm_l8sr,
  bandnm_l9sr];

//---------------------------------------------------------------------------//
// Step 0b: Define the direct inputs and run the full function
//---------------------------------------------------------------------------//
var y = 0; // refers to satellite sensor (count starts from 0). 0(S1), 2(S2SR), 4(L8SR), 5(L9SR)
            // to export S1 ARD images themselves, please see script 8gee_exports1ARD
var z = "all" //which points to get data for? 1 = first onboard group, 2 = second, "all" is all
var x = 0; //0=myanmar, 1=peru, 2=us

data(sat[y], filenames[y], allbands[y], z, disp[x], countryName[x]);


/*
//run the 'data' function for all satellites at once
var satindex = ee.List.sequence(0,5); //6 total sensors
satindex.size().evaluate(function(nlist){
  for(var y = 0; y<nlist; y++){
    data(country[x], sat[y], filenames[y], allbands[y], s1_despeckle[z]);
  }
});
*/
  
///////////////////////////////////////////////////////
///////////Other useful code ////////////////////////
///////////////////////////////////////////////////////

//1. Export each image as a separate csv
/*
//This loop exports a different chart for each of the elements in the list "loop".
//I don't need this now because I export everything together above, but 
//this is helpful to have.

// loop on client side
loop.size().evaluate(function(sizeCol){
  // loop on client side
  for (var i = 0; i<sizeCol;i++) {
    var toexport = ee.FeatureCollection(loop.get(i));
    var nam = "test_point_" + i;
   print('test '+i, toexport);
        // Export
  Export.table.toDrive(
    {collection: toexport,
    description: nam,
    folder: "dissertation",
    fileNamePrefix: nam,
    selectors: ['satellite', 'date', 'lat', 'lon', 'pointid', 
    'B1','B2','B3','B4','B5','B6','B7','B10','B11',
    'pixel_qa', 'radsat_qa', 'sr_aerosol']
    });
  };
});
*/

//2. Get separate chart for each image
/*
/// this uses test1 (imagecollection) from above
for (var i = 0; i < points.length; i++) {
  var coll = test1.get(i); 
  print(ui.Chart.image.series(coll, points, ee.Reducer.mean(), 30));
};
*/

/*
//3. function to mask clouds for landsat 8
function maskL8sr(image) {
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
  var cloudShadowBitMask = (1 << 3);
  var cloudsBitMask = (1 << 5);
  // Get the pixel QA band.
  var qa = image.select('pixel_qa');
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
               .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return image.updateMask(mask);
}
  
  
//4. function to get NDVI for landsat 8
var makeNDVI = function(image){
  var ndvi = image.normalizedDifference(['B5', 'B4']);
  return ndvi.set('system:time_start', image.get('system:time_start'));
};

*/

