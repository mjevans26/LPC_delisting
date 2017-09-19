/*
Generalized Change Detection Algorithms using Sentinel 2 data in Google Earth Engine
Version: 1.0
Date: 2017-09-12
License: CC-BY
Authors: Michael Evans (Defenders of Wildlife, Endangered Species Program)
*/
//DEFINE AREA AND TIME PERIODS
//Define area of interest (aoi)
//This is can be a polygon drawn on the fly in Earth Engine, or an uploaded shapefile
var aoi = geometry;

var today = ee.Date(Date.now());

var doi = ee.Date('2015-09-01');

var bands_re = ee.List(['B1', 'B2', 'B3', 'B4', 'B8', 'B11', 'B12', 'Pan', 'B10', 'TIR1', 'TIR2', 'BQA', 'fmask']);
var bands_orig = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12'];
var bands_class = ee.List(['B2', 'B3', 'B4', 'B8', 'B11', 'B12', 'ndvi', 'ndsi', 'nbr', 'system:time_start']);

//DEFINE FUNCTIONS FOR IMAGE PROCESSING
//Function to normalize band values across defined region
//Returns image with original bands normalized by global mean/sd values
// img: Single or multi-band image of reflectange values
// bnds: list of named bands to be normalized
// region: geometry constraining analysis to area of interest
var norm = function(img, bnds, region) {
  var imn = img.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: 30,
    maxPixels: 1000175326
  });
  var isd = img.reduceRegion({
    reducer: ee.Reducer.stdDev(),
    geometry: region,
    scale:30,
    maxPixels: 1000175326
  });

var array = img.select(bnds).toArray();
imn = imn.toArray(bnds);
isd = isd.toArray(bnds);

return array.subtract(imn).divide(isd).arrayProject([0]).arrayFlatten([bnds]);

};

//Function to calculate normalized difference metrics and add bands to collections
// img: Multi-band image of reflectange values
// NIR: band name corresponding to near infrared spectrum
var ND = function(img, NIR, R, G, SWIR1, SWIR2){
  var NBR = img.normalizedDifference([NIR, SWIR2]).rename(["nbr"]);
  var NDSI = img.normalizedDifference([G, SWIR1]).rename(["ndsi"]);
  //var NNDSI = norm(NDSI, ["ndsi"], region).rename(["nndsi"]);
  var NDVI = img.normalizedDifference([NIR, R]).rename(["ndvi"]);
  //var NNDVI = norm(NDVI, ["ndvi"], region).rename(["nndvi"]);

  return img.addBands(NDSI)
  //.addBands(NNDSI)
  .addBands(NDVI).addBands(NBR);
  //.addBands(NNDVI);
};

//Functions to mask clouds for S2 and LS8 data
// Convert processed sentinel toa reflectance to raw values,
//and extract azimuth/zenith metadata
//Returns multiband S2 image with values (0, 1]
// img: S2 MSI level 1-C image
function sentinel2toa(img) {
  var toa = img.select(['B1','B2','B3','B4', 'B5', 'B6','B7','B8', 'B8A','B9','B10', 'B11','B12'])  
     .divide(10000)
     .set('solar_azimuth',img.get('MEAN_SOLAR_AZIMUTH_ANGLE'))
     .set('solar_zenith', img.get('MEAN_SOLAR_ZENITH_ANGLE'));
  return toa;
}

//Function to identify clouds and shadows using
//classical Bayesian classifier from Hollstein et al. 2016
//top and bottom portion of distribution are shadows/clouds
//NOT USED CURRENTLY
var bayesP = function(toa){
var B1 = toa.select('B1');
var B2 = toa.select('B2');
var B3 = toa.select('B3');
var B8 = toa.select('B8A');
var B9 = toa.select('B9');
var B10 = toa.select('B10');
var B12 = toa.select('B12');
var Ib10b2 = B10.subtract(B2).divide(B10.add(B2));
var Ib2b8 = B2.subtract(B8).divide(B2.add(B8));
var index = B3.multiply(B9.subtract(B1)).multiply(B12).multiply(Ib10b2).multiply(Ib2b8);
return (index);
};

//Function to id clouds and cirrus using optimized classification tree
//from Hollstein et al. 2016
//Returns classified image with binary [0,1] 'cloud' and 'cirrus' bands
// toa: Sentinel 2 image converted to surface radiance
function ct(toa){
  var B1 = toa.select('B1');
  var B2 = toa.select('B2');
  var B3 = toa.select('B3');
  var B5 = toa.select('B5');
  var B6 = toa.select('B6');
  var B7 = toa.select('B7');
  var B8 = toa.select('B8A');
  var B9 = toa.select('B9');
  var B10 = toa.select('B10');
  var B11 = toa.select('B11');
  var B12 = toa.select('B12');
  var R15 = B1.divide(B5);
  var R29 = B2.divide(B9);
  var R210 = B2.divide(B10);
  var R511 = B5.divide(B11);
  var S1110 = B11.subtract(B10);
  var S67 = B6.subtract(B7);
  var S37 = B3.subtract(B7);
  var S911 = B9.subtract(B11);
  
  var cirrus1 = B3.lt(0.319).and(B8.gt(0.166)).and(R210.lt(14.689)).and(R29.gt(0.788));
  var cirrus2 = B3.gt(0.319).and(R511.lt(4.33)).and(S1110.lt(0.255)).and(S67.lt(0.016));
  var cirrus = cirrus1.or(cirrus2).rename(['cirrus']);
  
  var cloud1 = B3.gt(0.319).and(R511.lt(4.330)).and(S1110.lt(0.255)).and(S67.lt(0.016));
  var cloud2 = B3.gt(0.319).and(R511.lt(4.330)).and(S1110.gt(0.255)).and(B1.gt(0.3));
  var cloud = cloud1.or(cloud2).rename(['cloud']);
  
  //var shadow1 = B3.lt(0.319).and(B8.lt(0.166)).and(S37.lt(0.027)).and(S911.gt(0.097));
  //var shadow2 = B3.lt(0.319).and(B8.lt(0.166)).and(S37.gt(0.027)).and(S911.gt(0.021));
  //var shadow3 = B3.gt(0.319).and(R511.gt(4.330)).and(B3.lt(0.525)).and(R15.gt(1.184));
  //var shadow = shadow1.or(shadow2.or(shadow3)).rename(['shadow']);
 
 var shadow1 = B8.lt(0.181).and(B8.lt(0.051)).and(B9.lt(0.01)).and(B2.lt(0.073));
 var shadow2 = B8.lt(0.181).and(B8.lt(0.051)).and(B9.gt(0.01)).and(B3.gt(0.074));
 var shadow3 = B8.lt(0.181).and(B8.gt(0.051)).and(B12.lt(0.097)).and(B10.gt(0.011));
 var shadow4 = B8.lt(0.181).and(B8.gt(0.051)).and(B12.gt(0.097)).and(B10.gt(0.010));
 var shadow = shadow1.or(shadow2.or(shadow3.or(shadow4))).rename(['shadow']);
  return(ee.Image.cat(cirrus, cloud, shadow));
}

//Function to find dark pixels from threshold on sum of NIR, SWIR1, & SWIR2 bands
//Returns classified image with binary [0,1] 'dark' band
// toa: Sentinel 2 image converted to surface radiance
// thresh: threshold (0.2 - 0.5) value for sum of NIR, SWIR1 & SWIR2 bands
function findDarkPixels(toa, thresh) {
  var darkPixels = toa.select(['B8','B11','B12']).reduce(ee.Reducer.sum()).lt(thresh);
  var filtered = darkPixels.focal_mode(2, 'square', 'pixels');
  return filtered.rename(['dark']);
}

//Function to project cloud footprint based on sun azimuth and zentith
//Returns binary [0,1] image identifying areas of potential cloud shadow
// img: S2 radiance image with azimuth and zenith metadata
// cloudMask: image with binary [0,1] 'cloud' band
// cloudHeights: list of cloud altitudes over which to iterate
function project_clouds(img, cloudMask, cloudHeights){
  // get solar azimuth and zenith
  // var meanAzimuth  = img.get('MEAN_SOLAR_AZIMUTH_ANGLE');    // raw metadata
  // var meanZenith   = img.get('MEAN_SOLAR_ZENITH_ANGLE');     // raw metadata
  var meanAzimuth  = img.get('MEAN_SOLAR_AZIMUTH_ANGLE');  // recalculated for daily mosaics
  var meanZenith   = img.get('MEAN_SOLAR_ZENITH_ANGLE');   // recalculated for daily mosaics
  
  // convert to radians
  var azR  = ee.Number(meanAzimuth).multiply(Math.PI).divide(180.0).add(ee.Number(0.5).multiply(Math.PI));
  var zenR = ee.Number(0.5).multiply(Math.PI).subtract(ee.Number(meanZenith).multiply(Math.PI).divide(180.0));
  //var azR  = ee.Number(meanAzimuth).add(180).multiply(Math.PI).divide(180.0);
  //var zenR = ee.Number(meanZenith).multiply(Math.PI).divide(180.0);
  
  // get scale of image
  var nominalScale = cloudMask.projection().nominalScale();
  
  // find the shadows
  var shadows = cloudHeights.map(function(cloudHeight){
    cloudHeight = ee.Number(cloudHeight);
    var shadowCastedDistance = zenR.tan().multiply(cloudHeight);  // distance shadow is cast
    // var x = azR.cos().multiply(shadowCastedDistance).divide(nominalScale);//.round();  // x distance of shadow
    // var y = azR.sin().multiply(shadowCastedDistance).divide(nominalScale);//.round();  // y distance of shadow
    var x = azR.cos().multiply(shadowCastedDistance).divide(nominalScale);//.round();  // x distance of shadow
    var y = azR.sin().multiply(shadowCastedDistance).divide(nominalScale);//.round();  // y distance of shadow
    return cloudMask.changeProj(cloudMask.projection(), cloudMask.projection().translate(x, y));
  });
  var shadowMask = ee.ImageCollection.fromImages(shadows).max();
  
  // remove cloudMask from shadows
  shadowMask = shadowMask.and(cloudMask.not()).rename(['shadow']);
  return shadowMask;
}

//Function to find dark areas intersecting with cloud shadows.
//*Need to make more inclusive...Bayesian output?
//Returns image with binary [0,1], 'shadows' band
// darkPix: single band image indicating dark areas
// shadowPot: single band image identifying projected cloud footprings
function shadow_x_dark (darkPix, shadowPot){
  //var shadow = shadowPot.and(darkPix).rename('shadows');
  var darkareas = darkPix.updateMask(darkPix.eq(1)).reduceToVectors({
    geometry: geometry,
    maxPixels:1000000000
  });
  var shadows = shadowPot.updateMask(shadowPot.eq(1)).reduceToVectors({
    geometry: geometry,
    maxPixels: 1000000000
  });
  var intersect = darkareas.filterBounds(shadows);
  
  var shadow = intersect.reduceToImage(['label'], ee.Reducer.first());
  return intersect.rename(['shadows']);
}

//Function combining all steps of cloud/shadow masking.
//Can be mapped over image collection
//Returns image collection with masked cloud/shadow pixels
// img: S2 MSI level-1C image
function cloud_and_shadow_mask(img) {
  var toa = sentinel2toa(img);
  var cloud = ct(toa);
  var dark = findDarkPixels(toa, 0.35);
  var shadow = project_clouds(toa, cloud.select('cloud').eq(1), ee.List.sequence(100,500,100));
  var shadows = shadow_x_dark(dark, shadow);
  var mask = cloud.or(shadows);
  return img.updateMask(mask).addBands(img.metadata('system:time_start'));
}

//Alternative, simple function for masking clouds from S2 images
var S2maskClouds = function(img) {
  return img.updateMask(img.select('QA60')
  .lt(10)
  )
  .addBands(img.metadata('system:time_start')
  );
};

//Function to mask clouds from Landsat, using GEE algorithm
var LSmaskClouds = function(img) {
  var scored = ee.Algorithms.Landsat.simpleCloudScore(img);
  var masked = scored.updateMask(scored.select(['cloud']).lte(20));
  return masked.addBands(img.metadata('system:time_start'));
};

//Function to apply cloud mask, reduce image collection to single image,
//and clip to area of interest for S2 and LS8.
//Returns image composite
// img_col: Landsat or Sentinel 2 image collection
// region: geometry to which image is clipped
// type: c('median', 'recent') reducer to apply to image collection
var LSprocess = function(img_col, region, type){
  var area = img_col.filterBounds(region);
  var cloud = area.map(LSmaskClouds);
  var med;
  if (type == 'recent'){
    med = cloud.qualityMosaic('system:time_start')
    .clipToCollection(region);
  } else if (type == 'median'){
      med = cloud.median()
      .clipToCollection(region);
  }
  var rename = med;
  return ND(rename, 'B5', 'B4', 'B3', 'B6', 'B7', aoi);
  //var final = norm(nd, bands_class, aoi);
};

var S2process = function(img_col, region, type){
  var area = img_col.filterBounds(region);
  var cloud = area.map(cloud_and_shadow_mask);
  var med;
  if (type == 'recent'){
    med = cloud.qualityMosaic('system:time_start')
    .clipToCollection(region);
  } else if (type == 'median'){
      med = cloud.median()
      .clipToCollection(region);
  }
  var rename = med;
  return ND(rename, 'B8', 'B4', 'B3', 'B11', 'B12', aoi);
  //var final = norm(nd, bands_class, aoi);
};

//CREATE BEFORE AND AFTER IMAGES
var before_on = S2process(
  S2.filterDate('2016-05-01', today), aoi, 'median')
  .select(bands_class);
var after_off = S2process(
  S2.filterDate('2016-12-01', '2017-04-01'), aoi, 'median')
  .select(bands_class);
var before_off = S2process(
  S2.filterDate('2015-05-01', '2016-03-01'), aoi, 'median')
  .select(bands_class);
var after_on = S2process(
  S2.filterDate('2017-05-01', today), aoi, 'recent')
  .select(bands_class);

//Map.addLayer(before_on, {}, 'before');
//Map.addLayer(after_on, {}, 'after_on');
//Map.addLayer(before_off, {}, 'before_off');
//Map.addLayer(after_off, {}, 'after_off');

//CALCULATE CHANGE METRICS BETWEEN IMAGES
//Function calculating 'change vector'
//Returns image with 'cv' band
// b: before image
// a: after image
// bnds: list of band names (typically R, G, B, NIR, SWIR1, SWIR2)
var cv = function(b, a, bnds){
return(
  b.select(bnds)
                  .subtract(a
                            .select(bnds))
                  .pow(2)
.reduce(ee.Reducer.sum())
.rename(['cv'])
);  
};

//Function calculating 'relative change vector maximum'
//Returs image with 'rcvmax' band
// b: before image
// a: after image
// bnds: list of band names (typically R, G, B, NIR, SWIR1, SWIR2)
var rcvmax = function(b, a, bnds){
  return(
    b.select(bnds)
    .subtract(a.select(bnds))
    .divide(b.select(bnds)
     .max(a.select(bnds)))
    .pow(2)
.reduce(ee.Reducer.sum())
.rename(['rcvmax'])
    );
};

//Function calculating NDVI, NBR, & NDSI
//Returns image with 'ndvi', 'nbr', & 'ndsi' bands
// b: before image
// a: after image
// bnds: list of band names (typically R, G, B, NIR, SWIR1, SWIR2)
var d = function(b, a, bnds){
  return(
    b.select(bnds).subtract(a.select(bnds))
    );
};

//CREATE IMAGE WITH BANDS FOR CHANGE METRICS CV, RCV, NDVI, NBR, NDSI
//calculate cv from before and after images
var cv_on = cv(before_on, after_on, bands_orig);
var cv_off = cv(before_off, after_off, bands_orig);
//calculate rcv from before and after images
var rcv_on = rcvmax(before_on, after_on, bands_orig);
var rcv_off = rcvmax(before_off, after_off, bands_orig);
//calculate combined normalized difference metrics from before and after images
var diff_on = d(before_on, after_on, ['ndvi', 'nbr', 'ndsi']);
var diff_off = d(before_off, after_off, ['ndvi', 'nbr', 'ndsi']);
//combine cv, rcv, and normalized difference images into single image
var change_on = cv_on.addBands(rcv_on).addBands(diff_on);
var change = cv_off.addBands(rcv_off).addBands(diff_off);

//Map.addLayer(change_on, {}, "c_on");
//Map.addLayer(change, {}, "c");

//CALCULATE Z-SCORES FOR EACH METRIC
//Function to calclulate mean, min, max, sd of change values across change image.
//Returns dictionary with values for each band.
// img: spectral change image with bands ['cv', 'rcvmax', 'ndvi', 'ndsi', 'nbr']
// scl: scale (m) at which to sample image values.
var stats = function(img, scl){
  var mean = img.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: aoi,
  scale: scl,
  maxPixels: 100000000000
  }).rename(['cv', 'rcvmax', 'ndvi', 'ndsi', 'nbr'], ['cv_mn', 'rcv_mn', 'ndvi_mn', 'ndsi_mn', 'nbr_mn']);

  var sd = img.reduceRegion({
  reducer: ee.Reducer.stdDev(),
  geometry: aoi,
  scale: scl,
  maxPixels: 100000000000
  }).rename(['cv', 'rcvmax', 'ndvi', 'ndsi', 'nbr'], ['cv_sd', 'rcv_sd', 'ndvi_sd', 'ndsi_sd', 'nbr_sd']);

  var min = img.reduceRegion({
  reducer: ee.Reducer.min(),
  geometry: aoi,
  scale: scl,
  maxPixels: 100000000000
  }).rename(['cv', 'rcvmax', 'ndvi', 'ndsi', 'nbr'], ['cv_min', 'rcv_min', 'ndvi_min', 'ndsi_min', 'nbr_min']);  
  
  return mean.combine(sd).combine(min);
};

//Define function to select thresholds for change metrics.
//Returns binary [0,1] image selecting pixels exceeding all thresholds.
// img: spectral change image with bands ['cv', 'rcvmax', 'ndvi', 'ndsi', 'nbr']
// t1: z-score threshold for detecting change in cv metric
// t2: z-score threshold for detecting change in rcv metric
// t3: z-score threshold for detecting change in ndvi metric
// t4: z-score threshold for detecting change in nbr metric
// t5: z-score threshold for detecting change in ndsi metric
var thresh = function(img, t1, t2, t3, t4, t5){
  
  var dic = stats(img, 30);
 var cv_thresh = ee.Number(dic.get('cv_min'))
 .add(ee.Number(dic.get('cv_sd'))
 .multiply(t1));

 var rcv_thresh = ee.Number(dic.get('rcv_min'))
 .add(ee.Number(dic.get('rcv_sd'))
 .multiply(t2));

 var ndvi_thresh = ee.Number(dic.get('ndvi_mn'))
 .add(ee.Number(dic.get('ndvi_sd'))
 .multiply(t3));

 var nbr_thresh = ee.Number(dic.get('nbr_mn'))
 .add(ee.Number(dic.get('nbr_sd'))
 .multiply(t4));
 
var ndsi_thresh = ee.Number(dic.get('ndsi_mn'))
 .add(ee.Number(dic.get('ndsi_sd'))
 .multiply(t5));

 return(
   img.select(['cv']).gt(cv_thresh)
 .and(img.select(['rcvmax']).gt(rcv_thresh))
 .and(img.select(['ndvi']).gt(ndvi_thresh))
 //.and(img.select(['nbr']).gt(nbr_thresh))
 .and(img.select(['ndsi']).lt(ndsi_thresh))
);

};

//thresholds selecting values above critical ndsi
var thresh2 = function(img, t1, t2, t3, t4, t5){
  
  var dic = stats(img, 30);
 var cv_thresh = ee.Number(dic.get('cv_min'))
 .add(ee.Number(dic.get('cv_sd'))
 .multiply(t1));

 var rcv_thresh = ee.Number(dic.get('rcv_min'))
 .add(ee.Number(dic.get('rcv_sd'))
 .multiply(t2));

 var ndvi_thresh = ee.Number(dic.get('ndvi_mn'))
 .add(ee.Number(dic.get('ndvi_sd'))
 .multiply(t3));

 var nbr_thresh = ee.Number(dic.get('nbr_mn'))
 .add(ee.Number(dic.get('nbr_sd'))
 .multiply(t4));
 
var ndsi_thresh = ee.Number(dic.get('ndsi_mn'))
 .add(ee.Number(dic.get('ndsi_sd'))
 .multiply(t5));

 return(
   img.select(['cv']).gt(cv_thresh)
 .and(img.select(['rcvmax']).gt(rcv_thresh))
 .and(img.select(['ndvi']).gt(ndvi_thresh))
 //.and(img.select(['nbr']).gt(nbr_thresh))
 .and(img.select(['ndsi']).gt(ndsi_thresh))
);
};

//thresholds based on previous LDA analysis
var lda = function(img, int, ld1, ld2, ld3, ld4, ldt){
  
  var dic = stats(img, 100);
 var cv_z = img.select('cv').subtract(ee.Number(dic.get('cv_min')))
 .divide(ee.Number(dic.get('cv_sd')))
 .multiply(ld1);

 var rcv_z = img.select('rcvmax').subtract(ee.Number(dic.get('rcv_min')))
 .divide(ee.Number(dic.get('rcv_sd')))
 .multiply(ld2);

 var ndvi_z = img.select('ndvi').subtract(ee.Number(dic.get('ndvi_mn')))
 .divide(ee.Number(dic.get('ndvi_sd')))
 .multiply(ld3);

var ndsi_z = img.select('ndsi').subtract(ee.Number(dic.get('ndvi_mn')))
 .divide(ee.Number(dic.get('ndsi_sd')))
 .multiply(ld4);

var ld = cv_z.add(rcv_z).add(ndvi_z).add(ndsi_z).add(int);

 return(
   ld.gt(ldt)
);

};

//Calculate z-scores and select pixels based on thresholds to change metric image
var thresh_on = thresh(change_on, 2, 2, 1, 0, 1);
var thresh_off = thresh(change_off.clipToCollection(aoi), 2.5, 2.5, 0.20, 0, -0.5);

//Eliminate single pixels, and condense convoluted areas with majority filter
var opened = thresh_off.focal_mode({radius: 1, kernelType: 'square', units: 'pixels'})
.focal_mode({radius: 1, kernelType: 'square', units: 'pixels'});

//Convert smoothed change threshold raster to polygons
var thresh_poly = opened.eq(1).reduceToVectors({
  geometry: aoi,
  scale: 10,
  geometryType: 'polygon',
  eightConnected: true,
  //need to dynamically determine max pixels, or just set absurdly high
  maxPixels: 750000000000
});

//DEFINE SHAPE METRIC FUNCTIONS & ADD TO CHANGE POLYGONS
//solidity metric (0,1] (Jiao&Liu 2012)
var SOLI = function(ft){
  var soli = ft.area(5).divide(ft.convexHull().area(5));
  return soli;
};

//form factor (0,1] (Jiao&Liu 2012)
var FORM = function(ft){
  var form = ft.area(5).multiply(4).multiply(3.14).divide(ft.perimeter(5).pow(2));
  return form;
};

//shape index (0.88, Inf] (Jiao&Liu 2012)
var SI = function(ft){
  var si = ft.perimeter(5).divide(ft.area(5).sqrt().multiply(4));
  return si;
};

//elongation metric (Jiao&Liu 2012)
var ELON = function(ft){
 //create bounding rectangle, returns geometry
 var rect = ft.bounds().geometry();
 //get rectangle coordinates, returns 'Geo JSON style List' can't access indexed elements?
 var coords = rect.coordinates().get(0);
 //create multipoint geometry from coordinates
 var points = ee.Geometry.MultiPoint(coords);
 var point1 = ee.Geometry(points.geometries().get(0));
 var point2 = ee.Geometry(points.geometries().get(1));
 var point3 = ee.Geometry(points.geometries().get(2));
 var dist1 = point1.distance(point2);
 var dist2 = point2.distance(point3);
 var elon = (dist1.max(dist2)).divide(dist1.min(dist2));
 return elon;
};

//Function to calculate shape metrics and add as attributes to features.
//Should be mapped across feature collection returned by change detection
//Returns feature collection with added attributes
var metrics = function(ft){
 var form = FORM(ft);
 var soli = SOLI(ft);
 var elon = ELON(ft);
 var si = SI(ft);
 var perim = ft.perimeter(5);
 var area = ft.area(5);
 var ratio = perim.divide(area);
 var score = soli.add(form).add((si.subtract(0.88)).pow(-1));
 return ft.set({'form': form, 'soli': soli, 'elon': elon, 'si': si, 
                'perim': perim, 'area': area, 'ratio': ratio,
                'score': score
 });
};

//Map function calculating change metrics across all features in feature collection
var change_poly_met = thresh_poly.filter(ee.Filter.eq('label', 1)).map(metrics);

//Export image of threshold pixels
Export.image.toAsset({
  image: opened2,
  description: 'NewWells',
  scale: 10,
  maxPixels:750000000000
});

//Export polygons
Export.table.toDrive(change_poly_met, 'Well_Confirm_thresh');